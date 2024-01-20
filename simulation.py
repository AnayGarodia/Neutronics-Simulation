# 1. PREREQUISITES
# Importing required libraries
import openmc, openmc.deplete, os, plotly.graph_objects as go, plotly.io as pio, neutronics_material_maker as nmm, numpy as np, prettytable, math

# Constants
BLANKET_VARIATION_COUNT = int((500 - 100) / 25 + 1)
SOURCE_ENERGY_VARITATION_COUNT = int((14e6 - 0) / 1e6 + 1)
N_A = 6.02214076 * 10e23
dir_name = "/home/anay/conda/envs/paramak_env/fusion"
ALL_MATERIALS = nmm.AvailableMaterials()
BATCHES = 10
PARTICLES = 10000
NEUTRONS_PER_SECOND = 10e20
DATA = []

# Configuring cross sections
openmc.config["cross_sections"] = "/home/anay/sections/b/cross_sections.xml"
openmc.config["chain_file"] = "/home/anay/sections/b/chain.xml"

# 2. FUNCTIONS


# Removes all files with given extension from directory. File extension example = ".png"
def file_remover(file_extension: str, directory_name: str) -> None:
    directory = os.listdir(directory_name)
    for file in directory:
        if file.endswith(file_extension):
            os.remove(os.path.join(directory_name, file))


# Searches neutronic material maker library and returns a 2D Array containing materials
def create_materials(material_names, temperatures):
    materials_list = []
    for i in range(len(material_names)):
        temp = []
        for material in ALL_MATERIALS.keys():
            if not material.lower().find(material_names[i].lower()) == -1:
                temp.append(
                    nmm.Material.from_library(
                        name=material, temperature=temperatures[i]
                    ).openmc_material
                )
        materials_list.append(temp)

    return materials_list


def create_settings(particles: int, batches: int, source, inactive: int, run_mode: str):
    return openmc.Settings(
        particles=particles,
        batches=batches,
        source=source,
        inactive=inactive,
        run_mode=run_mode,
    )


def create_source(energy: int):
    return openmc.IndependentSource(
        space=openmc.stats.Point((0, 0, 0)),
        angle=openmc.stats.Isotropic(),
        energy=openmc.stats.Discrete([energy], [1]),
    )


# Returns geometry with neutron surface sphere radius being 1 m and thickness of both first wall and last wall being 2 cm.


def create_geometry(blanket_radius: int, materials):
    radii = [100, 102, blanket_radius + 102, blanket_radius + 104]

    global cells
    cells = []

    for i, radius in enumerate(radii):
        sphere = openmc.Sphere(r=radius)
        if i == len(radii) - 1:
            sphere.boundary_type = "vacuum"

        region = -sphere if i == 0 else +openmc.Sphere(r=radii[i - 1]) & -sphere

        cell = openmc.Cell(region=region)

        if i > 0:
            cell.fill = materials[i - 1]

        cells.append(cell)

    geometry = openmc.Geometry(cells)

    color = {
        cells[0]: "orange",
        cells[1]: "yellow",
        cells[2]: "blue",
        cells[3]: "green",
    }
    plot = geometry.plot(basis="xz", colors=color)
    file_remover("xz-cell.png", dir_name)
    plot.figure.savefig("xz-cell.png")

    return geometry


def create_filter(type: str, arg):
    if type == "cell":
        return openmc.CellFilter(arg)
    elif type == "particle":
        return openmc.ParticleFilter(arg)
    else:
        raise ValueError(
            "Invalid filter type. Type can be either 'Cell' or 'Particle'."
        )


def create_tally(name: str, filter, score):
    tal = openmc.Tally(name=name)
    tal.filters = filter
    tal.scores = score
    return tal


def sphere_vol(rad_out, rad_in):
    return 4 * math.pi / 3 * (math.pow(rad_out, 3) - math.pow(rad_in, 3))


# ASSUMPTIONS: TEMPERATURE FOR MATERIALS, 2*THRESHOLD ENERGY = 80.0, RECOMBINATION FACTOR = 0.8
def dpa_and_tbr():
    dpa_values = []
    tbr_values = []

    def calculate_dpa(damage_energy, material, volume):
        displacement_per_neutron_with_recombination = damage_energy * 0.8 / 80
        number_of_neutrons_per_year = NEUTRONS_PER_SECOND * 60 * 60 * 24 * 365.25
        displacement_for_all_atoms = (
            displacement_per_neutron_with_recombination * number_of_neutrons_per_year
        )
        number_of_atoms = material.density * volume * N_A / material.average_molar_mass
        return displacement_for_all_atoms / number_of_atoms

    mats = create_materials(
        ["steel, stainless 202", "lithium", "Steel, boron"], [1200, 1200, 1200]
    )
    for i in range(BLANKET_VARIATION_COUNT - 15):
        temp_dpa = []
        temp_tbr = []
        for j in range(len(mats[1]) - 9):
            blanket_thickness = 100 + i * 25
            materials = [mats[0][0], mats[1][j], mats[2][0]]
            geometry = create_geometry(blanket_thickness, materials)

            source = create_source(14.07e6)
            settings = create_settings(PARTICLES, BATCHES, source, 0, "fixed source")

            first_wall_filter = create_filter("cell", cells[1])
            blanket_filter = create_filter("cell", cells[2])
            last_wall_filter = create_filter("cell", cells[3])
            neutron_filter = create_filter("particle", ["neutron"])
            energy_filter = openmc.EnergyFilter.from_group_structure("CCFE-709")

            # DPA Tallies
            first_wall_tally = create_tally(
                "first_wall_tally", [first_wall_filter], ["444"]
            )
            last_wall_tally = create_tally(
                "last_wall_tally", [last_wall_filter], ["444"]
            )

            # TBR Tally
            blanket_tbr_tally = create_tally(
                "blanket_tbr_tally", [blanket_filter], ["(n,Xt)"]
            )
            blanket_tbr_tally.nuclides = ["Li6", "Li7"]

            first_wall_spectra_tally = create_tally(
                "first_wall_spectra_tally",
                [first_wall_filter, neutron_filter, energy_filter],
                ["flux"],
            )
            blanket_spectra_tally = create_tally(
                "blanket_spectra_tally",
                [blanket_filter, neutron_filter, energy_filter],
                ["flux"],
            )
            last_wall_spectra_tally = create_tally(
                "last_wall_spectra_tally",
                [last_wall_filter, neutron_filter, energy_filter],
                ["flux"],
            )

            tallies = openmc.Tallies(
                [
                    first_wall_tally,
                    last_wall_tally,
                    blanket_tbr_tally,
                    first_wall_spectra_tally,
                    blanket_spectra_tally,
                    last_wall_spectra_tally,
                ]
            )

            model = openmc.model.Model(geometry, materials, settings, tallies)

            file_remover(".h5", dir_name)

            results_filename = model.run()
            results = openmc.StatePoint(results_filename)

            first_wall_raw = results.get_tally(name="first_wall_tally")
            last_wall_raw = results.get_tally(name="last_wall_tally")
            tbr_tally = results.get_tally(name="blanket_tbr_tally")
            first_wall_flux_raw = results.get_tally(name="first_wall_spectra_tally")
            blanket_flux_raw = results.get_tally(name="blanket_spectra_tally")
            last_wall_flux_raw = results.get_tally(name="last_wall_spectra_tally")

            df1 = first_wall_raw.get_pandas_dataframe()
            df2 = last_wall_raw.get_pandas_dataframe()
            df3 = tbr_tally.get_pandas_dataframe()

            first_wall_damage_energy_in_ev = df1["mean"].sum()
            last_wall_damage_energy_in_ev = df2["mean"].sum()

            first_wall_volume = sphere_vol(102, 100)
            blanket_volume = sphere_vol(blanket_thickness + 102, 102)
            last_wall_volume = sphere_vol(
                blanket_thickness + 104, blanket_thickness + 102
            )

            first_wall_flux = (
                first_wall_flux_raw.mean.flatten()
                * NEUTRONS_PER_SECOND
                / first_wall_volume
            )
            blanket_flux = (
                blanket_flux_raw.mean.flatten() * NEUTRONS_PER_SECOND / blanket_volume
            )
            last_wall_flux = (
                last_wall_flux_raw.mean.flatten()
                * NEUTRONS_PER_SECOND
                / last_wall_volume
            )

            fig = go.Figure()

            # Add traces for each dataset
            fig.add_trace(
                go.Scatter(
                    x=energy_filter.values[:-1],
                    y=first_wall_flux,
                    mode="lines",
                    name="First Wall Flux",
                    line=dict(color="dodgerblue", width=1),
                )
            )

            fig.add_trace(
                go.Scatter(
                    x=energy_filter.values[:-1],
                    y=blanket_flux,
                    mode="lines",
                    name="Blanket Flux",
                    line=dict(color="seagreen", width=1),
                )
            )

            fig.add_trace(
                go.Scatter(
                    x=energy_filter.values[:-1],
                    y=last_wall_flux,
                    mode="lines",
                    name="Last Wall Flux",
                    line=dict(color="tomato", width=1),
                )
            )

            # Set the x and y axis to log scale
            fig.update_xaxes(type="log")
            fig.update_yaxes(type="log")

            # Set the x and y axis labels
            fig.update_layout(xaxis_title="Energy [eV]", yaxis_title="Flux [n/cm2-s]")

            # Show the figure
            pio.write_image(
                fig, f"png/FLUX_{materials[1].name}_{blanket_thickness}.png"
            )

            tbr_tally_result = df3["mean"].sum()

            # sums up all the values in the std. dev. column
            tbr_tally_std_dev = df3["std. dev."].sum()

            first_wall_dpa = calculate_dpa(
                first_wall_damage_energy_in_ev, materials[0], first_wall_volume
            )
            last_wall_dpa = calculate_dpa(
                last_wall_damage_energy_in_ev, materials[2], last_wall_volume
            )
            temp_dpa.append([first_wall_dpa, last_wall_dpa])
            temp_tbr.append(tbr_tally_result)

        dpa_values.append(temp_dpa)
        tbr_values.append(temp_tbr)

    blanket_sizes = list(
        range(100, 501, 25)
    )  # blanket sizes from 100 to 500 with a gap of 25

    materials = [
        material.name for material in mats[1]
    ]  # materials are members of mats[1]

    print(dpa_values)
    print(tbr_values)
    print(materials)

    # plot dpa values with a line plot using plotly y axis log x axis blanket size linear.
    for i in range(len(materials) - 9):
        tbr_vals = []
        dpa_vals1 = []
        dpa_vals2 = []
        for j in range(len(blanket_sizes) - 15):
            tbr_vals.append(tbr_values[j][i])
            dpa_vals1.append(dpa_values[j][i][0])
            dpa_vals2.append(dpa_values[j][i][1])

        fig_tbr = go.Figure()
        fig_dpa = go.Figure()

        fig_tbr.add_trace(
            go.Scatter(
                x=blanket_sizes,
                y=tbr_vals,
                mode="lines",
                name="TBR",
                line=dict(color="dodgerblue", width=1),
            )
        )

        # Set the x and y axis to log scale
        fig_tbr.update_xaxes(type="linear")
        fig_tbr.update_yaxes(type="linear")

        # Set the x and y axis labels
        fig_tbr.update_layout(
            xaxis_title="Blanket Thickness (in cm)", yaxis_title="TBR"
        )

        # Show the figure
        pio.write_image(fig_tbr, f"png/TBR_{materials[i]}.png")

        fig_dpa.add_trace(
            go.Scatter(
                x=blanket_sizes,
                y=dpa_vals1,
                mode="lines",
                name="DPA",
                line=dict(color="dodgerblue", width=1),
            )
        )

        fig_dpa.add_trace(
            go.Scatter(
                x=blanket_sizes,
                y=dpa_vals2,
                mode="lines",
                name="DPA",
                line=dict(color="tomato", width=1),
            )
        )

        # Set the x and y axis labels
        fig_dpa.update_layout(
            xaxis_title="Blanket Thickness (in cm)", yaxis_title="DPA"
        )

        # Show the figure
        pio.write_image(fig_dpa, f"png/DPA_{materials[i]}.png")

    DATA.append({"DPA": dpa_values})
    DATA.append({"TBR": tbr_values})


def depletion():
    mats = create_materials(
        ["steel, stainless 202", "lithium-lead", "Steel, boron"], [300, 300, 300]
    )
    mats[0][0].depletable = True
    mats[1][0].depletable = True
    mats[2][0].depletable = True

    materials = [mats[0][0], mats[1][0], mats[2][0]]
    for i in range(BLANKET_VARIATION_COUNT):
        blanket_thickness = 100 + i * 25
        geometry = create_geometry(blanket_thickness, materials)
        source = create_source(14.07e6)
        source.particles = "neutron"
        settings = create_settings(500, 2, source, 0, "fixed source")
        model = openmc.model.Model(geometry, materials, settings)
        flux_in_each_group, micro_xs = openmc.deplete.get_microxs_and_flux(
            model=model,
            domains=[cells[2]],
            energies="CCFE-709",
        )


if __name__ == "__main__":
    file_remover(".png", dir_name)
    file_remover(".png", f"{dir_name}/png")
    dpa_and_tbr()
