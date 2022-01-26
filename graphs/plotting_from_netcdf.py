#!/usr/bin/env python3

import argparse
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import itertools
import re


def plot_flux_surface(df: xr.Dataset, inputs: xr.Dataset, step=None):
    if step is None:
        step = df["rpts"].shape[1] // 20

    grid_spec = plt.GridSpec(2, 2, height_ratios=(1, 9))
    fig = plt.figure(figsize=(11.69, 8.27))
    title_ax = fig.add_subplot(grid_spec[0, :])
    ax0 = fig.add_subplot(grid_spec[1, 0])
    ax1 = fig.add_subplot(grid_spec[1, 1])

    title_ax.text(
        0.5,
        0.75,
        f"{df.software_name} {df.software_version}",
        fontsize="xx-large",
        horizontalalignment="center",
    )
    title_ax.text(
        0.5,
        0.25,
        f"Created:{df.date_created}\nRun ID: {df.id}",
        horizontalalignment="center",
    )
    title_ax.axis("off")

    ax0.plot(df["rpts"][:, ::step], df["zpts"][:, ::step], color="black", linewidth=1)

    ax0.axis("equal")

    ax0.set_ylabel("Z (m)")
    ax0.set_xlabel("R (m)")
    ax0.set_title("Flux surfaces")

    ax1.set_title("Input parameters")
    ax1.axis("off")

    all_inputs = ", ".join(
        [f"{key} = {value.data:.5G}" for key, value in inputs.data_vars.items()]
    )
    columned_inputs = re.sub(
        ", ",
        lambda m, c=itertools.count(1): m.group() if next(c) % 3 else "\n",
        all_inputs,
    )
    ax1.text(0.05, 0.5, columned_inputs, verticalalignment="center")

    fig.tight_layout()
    return fig, (title_ax, ax0, ax1)


def profiles(df: xr.Dataset):

    fig, ax = plt.subplots(ncols=2, nrows=2)

    ax[0, 0].plot(df.R, df.Ne, "r", label="$n_e$", linewidth=2)
    ax[0, 0].plot(df.R, df.Ni, "b", label="$n_i$", linewidth=2)
    ax[0, 0].legend(loc=1)
    ax[0, 0].xlabel("R (m)", fontsize=20)
    ax[0, 0].ylabel("Density $(m^{-3})$", fontsize=20)
    ax[0, 0].title("Electron density", fontsize=24)
    ax[0, 0].xticks(size=20)
    ax[0, 0].yticks(size=20)

    ax[0, 1].plot(df.R, df.B_t, label=r"$B_\phi$")
    ax[0, 1].plot(df.R, df.B_p, label=r"$B_\theta$")
    ax[0, 1].plot(df.R, df.B, label="$B$")
    ax[0, 1].xlabel("R (m)")
    ax[0, 1].ylabel("Field (T)")
    ax[0, 1].title("Magnetic fields")
    ax[0, 1].legend(loc=1)

    ax[1, 0].plot(df.R, df.q)
    ax[1, 0].title("Safety Factor")
    ax[1, 0].xlabel("R (m) ")
    ax[1, 0].ylabel("Safety factor")

    ax[0, 1].plot(df.R, df.T, "r--", label="$Temp$", linewidth=2)
    ax[0, 1].set_xlabel("R (m)", fontsize=20)
    ax[0, 1].set_ylabel("Temp (keV)", fontsize=20)
    ax[0, 1].tick_params("y", colors="r")
    ax[0, 1].xticks(size=20)
    ax[0, 1].yticks(size=20)
    ax2 = ax[0, 1].twinx()
    ax2.plot(df.R, df.p, "b.", label="$Pressure$")
    ax2.set_ylabel("Pressure", color="b")
    ax2.tick_params("y", colors="b")
    ax[0, 1].title("Electron Temperature", fontsize=24)

    return fig, ax


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        "plot_scene", description="Plot various quantities from SCENE"
    )
    parser.add_argument("filename", type=str, help="Filename of input netCDF file")
    parser.add_argument(
        "-o", "--output", type=str, help="Filename to save graphs to", default=None
    )

    args = parser.parse_args()

    df = xr.open_dataset(args.filename)
    inputs = xr.open_dataset(args.filename, group="inputs")

    fig_fluxsurface, _ = plot_flux_surface(df, inputs)

    plt.show()
    # fig_profiles, _ = profiles(df)

    if args.output is not None:
        with PdfPages(args.output) as pdf:
            pdf.savefig(fig_fluxsurface)
            # pdf.savefig(fig_profiles)
