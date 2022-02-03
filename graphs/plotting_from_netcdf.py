#!/usr/bin/env python3

import argparse
import xarray as xr
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from typing import List
import textwrap


A4_SIZE = (11.69, 8.27)


def plot_scene_title(title_ax: plt.Axes, df: xr.Dataset):
    title_ax.text(
        0.5,
        0.75,
        f"{df.title}",
        fontsize="xx-large",
        horizontalalignment="center",
        verticalalignment="bottom",
    )
    title_ax.text(
        0.5,
        0.5,
        f"{df.software_name} {df.software_version}",
        horizontalalignment="center",
        verticalalignment="center",
    )
    title_ax.text(
        0.5,
        0.0,
        f"Created: {df.date_created}\nRun ID: {df.id}",
        horizontalalignment="center",
        verticalalignment="top",
    )
    title_ax.axis("off")
    return title_ax


def plot_flux_surface(
    df: xr.Dataset, inputs: xr.Dataset, impurities: List[xr.Dataset], step=None
):
    if step is None:
        step = df["rpts"].shape[1] // 20

    grid_spec = plt.GridSpec(2, 3, height_ratios=(1, 9))
    fig = plt.figure(figsize=A4_SIZE)
    title_ax = fig.add_subplot(grid_spec[0, :])
    ax0 = fig.add_subplot(grid_spec[1:, 0])
    ax1 = fig.add_subplot(grid_spec[1, 1])
    ax2 = fig.add_subplot(grid_spec[1, 2])

    plot_scene_title(title_ax, df)

    ax0.plot(df["rpts"][:, ::step], df["zpts"][:, ::step], color="black", linewidth=1)
    ax0.plot(df["rpts"][:, 0], df["zpts"][:, 0], color="black", linewidth=4)

    ax0.axis("equal")

    ax0.set_ylabel("Z (m)")
    ax0.set_xlabel("R (m)")
    ax0.set_title("Flux surfaces")

    ax1.set_title("Input parameters")
    ax1.axis("off")

    all_inputs = "  ".join(
        [f"{key}\xa0=\xa0{value.data:.5G}" for key, value in inputs.data_vars.items()]
    )
    columned_inputs = textwrap.fill(all_inputs, width=35)
    ax1.text(0.05, 1, columned_inputs, verticalalignment="top")

    ax2.set_title("Impurity information")
    ax2.axis("off")

    def format_impurity(impurity) -> str:
        def fmt(attr: str) -> str:
            return f"{impurity[attr].data:.5G}"

        def fmt_key(attr: str) -> str:
            return f"{attr} = {fmt(attr)}"

        text = []
        text.append(f"Z = {fmt('Z')} M/Mp = {fmt('M')}")
        text.append(
            "  ".join([fmt_key(attr) for attr in ["at", "T0", "Ta", "Tped", "Tedg"]])
        )
        text.append(
            "  ".join([fmt_key(attr) for attr in ["an", "n0", "na", "nped", "nedg"]])
        )
        return "\n".join(text)

    columned_impurities = "\n\n".join(
        [format_impurity(impurity) for impurity in impurities]
    )
    ax2.text(
        0.05, 1, f"Number of impurities = {inputs.nimp.data}", verticalalignment="top"
    )
    ax2.text(0.05, 0.95, columned_impurities, verticalalignment="top")

    fig.tight_layout()
    return fig, (title_ax, ax0, ax1)


def profiles(df: xr.Dataset, min_R: float, max_R: float):

    fig, ax = plt.subplots(ncols=2, nrows=2, figsize=A4_SIZE)

    density_ax, field_ax, safety_ax, temperature_ax = ax.flat

    colour_cycle = [c for c in plt.rcParams["axes.prop_cycle"]]
    first_colour = colour_cycle[0]["color"]
    second_colour = colour_cycle[1]["color"]

    density_ax.plot(df.R, df.n_e, label="$n_e$", linewidth=2)
    density_ax.plot(df.R, df.n_i, label="$n_i$", linewidth=2)
    density_ax.legend(loc=1)
    density_ax.set_xlabel("R [m]")
    density_ax.set_ylabel("Density $[m^{-3}]$")
    density_ax.set_title("Electron density")

    field_ax.plot(df.R, df.B_T, label=r"$B_\phi$", linewidth=2)
    field_ax.plot(df.R, df.B_P, label=r"$B_\theta$", linewidth=2)
    field_ax.plot(df.R, df.B, label="$B$")
    field_ax.set_xlabel("R [m]")
    field_ax.set_ylabel("Field [T]")
    field_ax.set_title("Magnetic fields")
    field_ax.legend(loc=1)

    safety_ax.plot(df.R, df.q, label="q", linewidth=2, color=first_colour)
    safety_ax.set_title("Safety Factor")
    safety_ax.set_xlabel("R [m]")
    safety_ax.set_ylabel("Safety factor", color=first_colour)
    safety_ax.tick_params("y", colors=first_colour)

    collisions_ax = safety_ax.twinx()
    collisions_ax.plot(
        df.R, df.nu_star_e, "--", label=r"$\nu_{*e}$", linewidth=2, color=second_colour
    )
    collisions_ax.set_ylabel("Collisionality", color=second_colour)
    collisions_ax.tick_params("y", colors=second_colour)

    lines = safety_ax.get_lines() + collisions_ax.get_lines()
    labels = [line.get_label() for line in lines]
    collisions_ax.legend(lines, labels)

    temperature_ax.plot(
        df.R, df.T_e, "--", label="$T_e$", linewidth=2, color=first_colour
    )
    temperature_ax.plot(
        df.R, df.T_i, ".", label="$T_i$", linewidth=2, color=first_colour
    )
    temperature_ax.set_xlabel("R [m]")
    temperature_ax.set_ylabel("Temperature [keV]", color=first_colour)
    temperature_ax.tick_params("y", colors=first_colour)
    temperature_ax.set_title("Electron temperature and pressure")

    pressure_ax = temperature_ax.twinx()
    pressure_ax.plot(df.R, df.P_e, "-.", label="Pressure", color=second_colour)
    pressure_ax.set_ylabel("Pressure", color=second_colour)
    pressure_ax.tick_params("y", colors=second_colour)

    lines = temperature_ax.get_lines() + pressure_ax.get_lines()
    labels = [line.get_label() for line in lines]
    pressure_ax.legend(lines, labels)

    for axis in list(ax.flat) + [pressure_ax]:
        axis.set_xlim(min_R, max_R)
        axis.set_ylim(bottom=0)

    fig.tight_layout()

    return fig, (ax, pressure_ax)


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

    impurities = [
        xr.open_dataset(args.filename, group=f"inputs/impurity_{i:02d}")
        for i in range(1, inputs.nimp.data + 1)
    ]

    r_profiles = xr.open_dataset(args.filename, group="R_profiles")
    r_profiles["R"] = df.R

    fig_fluxsurface, _ = plot_flux_surface(df, inputs, impurities)

    fig_profiles, _ = profiles(r_profiles, df.rpts[:, :-1].min(), df.rpts.max())
    plt.show()

    if args.output is not None:
        with PdfPages(args.output) as pdf:
            pdf.savefig(fig_fluxsurface)
            pdf.savefig(fig_profiles)
