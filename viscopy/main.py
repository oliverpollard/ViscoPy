import numpy as np
import pandas as pd
from pathlib import Path
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt


MODEL_COLS = ["rad", "rhod", "vel_p", "vel_s", "gacc", "visc"]
DEFAULT_MODEL_FILE = Path(__file__).parent / "../resources/model.prem.l80.ump5.lm5"
CORE_LM_BOUNDARY = 3485500
LM_UM_BOUNDARY = 5700000
LITH_VISC = 1e43
EARTH_RADIUS = 6371000

MAXWELL_SOURCE = Path(__file__).parent / "../maxwell/maxwell_col_512.f"
MAXWELL_MODEL = Path(__file__).parent / "../maxwell/maxwell_col_512.x"
MAXWELL_SCRIPT = Path(__file__).parent / "../maxwell/run_maxwell.sh"


def plot_model(df, ax=None, xlim=None, type="line"):
    rads = df["rad"].values
    dist = (rads[-1] - rads) / 1000
    discont = list(np.diff(dist) == 0) + [False]

    def plot(ax, x, y, type):
        if type == "line":
            ax.plot(x, y)
        elif type == "scatter":
            ax.scatter(x, y)

    if ax is None:
        fig, ax = plt.subplots(2, 3, figsize=[20, 10])
        ax = ax.flatten()
        ax = {MODEL_COLS[idx]: ax[idx] for idx in range(6)}

    for idx, item in enumerate(df):
        if item == "visc":
            plot(ax[item], dist, np.log10(df[item]), type)
        else:
            plot(ax[item], dist, df[item], type)
        for d in dist[discont]:
            ax[item].axvline(d, ls="--", color="grey")
        ax[item].set_title(item)
        if xlim is None:
            xlim = [dist.max(), 0]
        ax[item].set_xlim(xlim)

    return ax


def get_boundary_value_positions(array_x, array_y, boundary):

    array_x, array_y = np.asarray(array_x), np.asarray(array_y)

    # outside of bounds: less than
    if boundary <= array_x.min():
        left_idx, right_idx = array_x.argmin(), array_x.argmin()
    # outside of bounds: greater than
    elif boundary >= array_x.max():
        left_idx, right_idx = array_x.argmax(), array_x.argmax()
    # inside bounds
    else:
        where_search = np.where(array_x == boundary)[0]
        # inside bounds: value already in array, and is a boundary
        if len(where_search) == 2:
            left_idx, right_idx = where_search[0], where_search[1]

        # inside bounds: value already in array, but is not a boundary
        elif len(where_search) == 1:
            left_idx, right_idx = where_search[0], where_search[0]

        # inside bounds: value not in array, get values either side
        else:
            # closest value to boundary
            min_idx = np.abs(array_x - boundary).argmin()
            closest_x = array_x[min_idx]

            if closest_x < boundary:
                left_x = closest_x
                # check if this is a boundary
                search_idx = np.where(array_x == left_x)[0]
                if len(search_idx) == 2:
                    # this is a boundary, and since it's too the left, take the most right value
                    left_idx = search_idx[-1]
                else:
                    left_idx = search_idx[0]
                right_idx = left_idx + 1
                right_x = array_x[right_idx]

            else:
                right_x = closest_x
                # check if this is a boundary
                search_idx = np.where(array_x == right_x)[0]
                if len(search_idx) == 2:
                    # this is a boundary, and since it's too the right, take the most left value
                    right_idx = search_idx[0]
                else:
                    right_idx = search_idx[0]
                left_idx = right_idx - 1
                left_x = array_x[left_idx]

    return left_idx, right_idx


def generate_visc_model(visc_lm, visc_um, lith_thickness, save=None):
    um_lith_boundary = EARTH_RADIUS - lith_thickness

    visc_model = {
        visc_lm: [CORE_LM_BOUNDARY, LM_UM_BOUNDARY],
        visc_um: [LM_UM_BOUNDARY, um_lith_boundary],
        LITH_VISC: [um_lith_boundary, EARTH_RADIUS],
    }

    return visc_model


def generate_r_dist(r_dist_model=None):
    if r_dist_model is None:
        r_dist_model = [
            [(0, LM_UM_BOUNDARY), 35 * 1000],
            [(LM_UM_BOUNDARY, 6350000.0), 25 * 1000],
            [(6350000.0, EARTH_RADIUS), 10500],
        ]

    r_dist = []
    for (lower_r_dist, upper_r_dist), step in r_dist_model:
        r_dist.append(np.arange(lower_r_dist, upper_r_dist, step))
    r_dist.append(np.array([EARTH_RADIUS]))
    r_dist = np.unique(np.concatenate(r_dist))

    return r_dist


def load_model(model_path):
    model_data_array = np.genfromtxt(model_path)
    # radius, density, p-wave vel, s-wave vel, grav accel, viscosity
    model_df = pd.DataFrame(data=model_data_array, columns=MODEL_COLS)
    return model_df


def save_model(df, save):
    np.savetxt(save, df[MODEL_COLS].values, fmt=["%16.7e"] + ["%15.7e"] * 5)


class ModelGenerator:
    def __init__(self, model_df=None):
        if model_df is None:
            model_df = load_model(DEFAULT_MODEL_FILE)
        self.model_df = model_df

    @property
    def r_dist(self):
        return self.model_df["rad"].values

    @property
    def boundaries(self):
        discont = list(np.diff(self.r_dist) == 0) + [False]
        boundaries = list(self.r_dist[discont])
        return boundaries

    @property
    def r_dist_max(self):
        r_dist_max = self.r_dist.max()
        return r_dist_max

    def generate(self, visc_model, r_dist_new=None):
        if r_dist_new is None:
            r_dist_new = generate_r_dist()

        data_new = {}
        for item in ["rhod", "vel_p", "vel_s", "gacc", "visc"]:
            array_y = self.model_df[item].values

            f = interp1d(self.r_dist, array_y)

            r_dist_new_temp = r_dist_new.copy()
            array_y_new = f(r_dist_new_temp)

            for boundary in self.boundaries:
                # find position of boundary in original data, and get values
                boundary_vals = array_y[np.where(self.r_dist == boundary)[0]]

                # if the boundary is already in x_new, remove it to replace values with boundary values
                boundary_new_r = np.where(r_dist_new_temp == boundary)[0]
                if len(boundary_new_r) != 0:
                    r_dist_new_temp, array_y_new = np.delete(
                        r_dist_new_temp, boundary_new_r
                    ), np.delete(array_y_new, boundary_new_r)

                # find where to insert values in new array
                ii = np.searchsorted(r_dist_new_temp, boundary)
                r_dist_new_temp = np.insert(r_dist_new_temp, ii, [boundary, boundary])
                array_y_new = np.insert(array_y_new, ii, boundary_vals)

            data_new[item] = array_y_new

        r_dist_new = r_dist_new_temp

        visc_boundaries = []
        for _, visc_bounds in visc_model.items():
            for boundary in visc_bounds:
                if (
                    (boundary not in visc_boundaries)
                    and (boundary not in self.boundaries)
                    and (boundary != self.r_dist_max)
                ):
                    visc_boundaries.append(boundary)

        for item in ["rhod", "vel_p", "vel_s", "gacc", "visc"]:
            r_dist_new_temp = r_dist_new.copy()
            array_y_temp = data_new[item].copy()
            for boundary in visc_boundaries:
                (
                    left_boundary_val_idx,
                    right_boundary_val_idx,
                ) = get_boundary_value_positions(
                    r_dist_new_temp, array_y_temp, boundary
                )
                boundary_vals = (
                    array_y_temp[left_boundary_val_idx],
                    array_y_temp[right_boundary_val_idx],
                )
                ii = np.searchsorted(r_dist_new_temp, boundary)
                array_y_temp = np.insert(array_y_temp, ii, boundary_vals)
                r_dist_new_temp = np.insert(r_dist_new_temp, ii, [boundary, boundary])

            data_new[item] = array_y_temp

        r_dist_new = r_dist_new_temp

        for value, visc_bounds in visc_model.items():
            mod_idx = np.where(
                (r_dist_new >= visc_bounds[0]) & (r_dist_new <= visc_bounds[1])
            )[0][1:-1]
            data_new["visc"][mod_idx] = value

        data_new["rad"] = r_dist_new
        df_new = pd.DataFrame(data_new)

        return df_new

    @classmethod
    def from_model(cls, model_file):
        model_df = load_model(model_file)
        return cls(model_df=model_df)


def gen_maxwell_cmd(input_file, output_file, run_file=None):
    cmd = f"bash {MAXWELL_SCRIPT} {MAXWELL_MODEL} {str(input_file)} {str(output_file)}"

    if run_file is not None:
        with open(run_file, "w") as fp:
            fp.write(f"{cmd}\n")
    else:
        return cmd


def output_earth_model(visc_lm, visc_um, lith_thickness, output_file):
    mg = ModelGenerator()
    visc_model = generate_visc_model(
        visc_lm=visc_lm,
        visc_um=visc_um,
        lith_thickness=lith_thickness,
    )
    model = mg.generate(visc_model=visc_model)
    save_model(model, save=output_file)


def output_earth_models(sample_df, output_dir):
    output_dir = Path(output_dir)
    for idx in range(len(sample_df)):
        output_earth_model(
            visc_lm=sample_df.iloc[idx]["visc_lm"],
            visc_um=sample_df.iloc[idx]["visc_um"],
            lith_thickness=sample_df.iloc[idx]["lith_thickness"],
            output_file=output_dir / f"{idx}.emtxt",
        )


def gen_maxwell_cmds(input_files, run_file=None):
    cmds = []
    for file in input_files:
        input_file = Path(file).resolve()
        output_file = input_file.with_suffix(".em")
        cmd = gen_maxwell_cmd(
            input_file=input_file, output_file=output_file, run_file=None
        )
        cmds.append(cmd)

    if run_file is not None:
        with open(run_file, "w") as fp:
            for cmd in cmds:
                fp.write(f"{cmd}\n")
    else:
        return cmds
