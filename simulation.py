# 1. PREREQUISITES
# Importing required libraries
import openmc, openmc.deplete, os, plotly.graph_objects as go, plotly.io as pio, neutronics_material_maker as nmm, math, pandas as pd

# Constants
BLANKET_VARIATION_COUNT = int((150 - 10) / 5 + 1)
N_A = 6.02214076 * 10e23
dir_name = "/home/anay/conda/envs/paramak_env/fusion"
ALL_MATERIALS = nmm.AvailableMaterials()
BATCHES = 10
PARTICLES = 25000
NEUTRONS_PER_SECOND = 10e20

# Configuring cross sections
openmc.config["cross_sections"] = "/home/anay/sections/b/cross_sections.xml"


# Removes all files with given extension from directory. File extension example = ".png"
def file_remover(file_extension: str, directory_name: str) -> None:
    directory = os.listdir(directory_name)
    for file in directory:
        if file.endswith(file_extension):
            os.remove(os.path.join(directory_name, file))


def create_materials(material_names, temperatures):
    materials_list = []
    for i in range(len(material_names)):
        temp = []
        if material_names[i] == "lithium":
            matt = nmm.Material.from_library(
                name="Lithium",
                temperature=temperatures[i],
                enrichment=40,
                enrichment_target="Li6",
                enrichment_type="ao",
            ).openmc_material
            matt.name = "Enriched Lithium"
            matt1 = nmm.Material.from_library(
                name="lithium-lead",
                temperature=temperatures[i],
                enrichment=40,
                enrichment_target="Li6",
                enrichment_type="ao",
            ).openmc_material
            matt1.name = "Enriched Lead Lithium"
            matt2 = nmm.Material.from_library(
                name="lithium-lead", temperature=temperatures[i]
            ).openmc_material
            matt2.name = "Lead Lithium"
            temp = [
                matt2,
                nmm.Material.from_library(
                    name="Lithium", temperature=temperatures[i]
                ).openmc_material,
                matt,
                matt1,
                nmm.Material.from_library(
                    name="Li8PbO6", temperature=temperatures[i]
                ).openmc_material,
                nmm.Material.from_library(
                    name="FLiBe", temperature=temperatures[i]
                ).openmc_material,
                nmm.Material.from_library(
                    name="FLiNaBe", temperature=temperatures[i]
                ).openmc_material,
                nmm.Material.from_library(
                    name="Li4SiO4", temperature=temperatures[i]
                ).openmc_material,
                nmm.Material.from_library(
                    name="Li2ZrO3", temperature=temperatures[i]
                ).openmc_material,
                nmm.Material.from_library(
                    name="Li2TiO3", temperature=temperatures[i]
                ).openmc_material,
            ]
        else:
            for material in ALL_MATERIALS.keys():
                if not material.lower().find(material_names[i].lower()) == -1:
                    temp.append(
                        nmm.Material.from_library(
                            name=material, temperature=temperatures[i], pressure=101325
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


def create_geometry(blanket_radius: int, materials, n: int):
    radii = [100, 102]
    if n == 4:
        radii.append(blanket_radius + 102)
        radii.append(blanket_radius + 104)
    if n == 6:
        radii.append(107)
        radii.append(blanket_radius + 107)
        radii.append(blanket_radius + 112)
        radii.append(blanket_radius + 114)
    if n == 7:
        radii.append(107)
        radii.append(blanket_radius + 107)
        radii.append(blanket_radius + 112)
        radii.append(blanket_radius + 117)
        radii.append(blanket_radius + 119)

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

    if blanket_radius == 50:

        if n == 4:
            color = {
                cells[0]: "lightgray",  # Dark Slate Gray
                cells[1]: "royalblue",  # Royal Blue
                cells[2]: "seagreen",  # Tomato
                cells[3]: "black",  # Sea Green
            }

        if n == 6:
            color = {
                cells[0]: "lightgray",  # Dark Slate Gray
                cells[1]: "royalblue",  # Royal Blue
                cells[2]: "tomato",  # Tomato
                cells[3]: "seagreen",  # Sea Green
                cells[4]: "indigo",  # Indigo
                cells[5]: "black",  # Saddle Brown
            }

        if n == 7:
            color = {
                cells[0]: "lightgray",  # Dark Slate Gray
                cells[1]: "royalblue",  # Royal Blue
                cells[2]: "tomato",  # Tomato
                cells[3]: "seagreen",  # Sea Green
                cells[4]: "indigo",  # Indigo
                cells[5]: "saddlebrown",  # Saddle Brown
                cells[6]: "black",  # Dark Orchid
            }

        plot = geometry.plot(basis="xz", colors=color)
        plot.figure.savefig(f"xz-cell{n}.png")

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


def dpa_and_tbr():
    def calculate_dpa(damage_energy, material, volume):
        displacement_per_neutron_with_recombination = damage_energy * 0.8 / 80
        number_of_neutrons_per_year = NEUTRONS_PER_SECOND * 60 * 60 * 24 * 365.25
        displacement_for_all_atoms = (
            displacement_per_neutron_with_recombination * number_of_neutrons_per_year
        )
        number_of_atoms = material.density * volume * N_A / material.average_molar_mass
        return displacement_for_all_atoms / number_of_atoms

    for q in range(3, 4):
        dpa_values = []
        tbr_values = []
        leakage_fraction = []
        mats = []
        if q == 0:
            mats = create_materials(
                [
                    "steel, stainless 202",
                    "lithium",
                    "Steel, boron",
                ],
                [1200, 1200, 1200],
            )
        if q == 1:
            mats = create_materials(
                [
                    "steel, stainless 202",
                    "Pb",
                    "lithium",
                    "Titanium Hydride",
                    "Steel, boron",
                ],
                [1200, 1200, 1200, 1200, 1200],
            )
        if q == 2:
            mats = create_materials(
                [
                    "steel, stainless 202",
                    "Pb",
                    "lithium",
                    "H2O",
                    "Steel, boron",
                ],
                [1200, 1200, 1200, 1200, 1200],
            )
        if q == 3:
            mats = create_materials(
                [
                    "steel, stainless 202",
                    "Pb",
                    "lithium",
                    "H2O",
                    "Titanium Hydride",
                    "Steel, boron",
                ],
                [1200, 1200, 1200, 1200, 1200, 1200],
            )

        for i in range(BLANKET_VARIATION_COUNT):
            temp_dpa = []
            temp_tbr = []
            temp_leakage = []
            x = 2
            if q == 0:
                x = 1
            for j in range(len(mats[x])):
                blanket_thickness = 10 + i * 5
                materials = []
                if q == 0:
                    materials = [
                        mats[0][0],
                        mats[1][j],
                        mats[2][0],
                    ]
                if q == 1 or q == 2:
                    materials = [
                        mats[0][0],
                        mats[1][0],
                        mats[2][j],
                        mats[3][0],
                        mats[4][0],
                    ]
                if q == 3:
                    materials = [
                        mats[0][0],
                        mats[1][0],
                        mats[2][j],
                        mats[3][0],
                        mats[4][0],
                        mats[5][0],
                    ]
                if q == 0:
                    x = 4
                elif q == 1 or q == 2:
                    x = 6
                elif q == 2:
                    x = 6
                elif q == 3:
                    x = 7
                geometry = create_geometry(blanket_thickness, materials, x)

                source = create_source(14.07e6)
                settings = create_settings(
                    PARTICLES, BATCHES, source, 0, "fixed source"
                )

                if q == 0:
                    first_wall_filter = create_filter("cell", cells[1])
                    blanket_filter = create_filter("cell", cells[2])
                    last_wall_filter = create_filter("cell", cells[3])
                elif q == 1:
                    first_wall_filter = create_filter("cell", cells[1])
                    multiplier_filter = create_filter("cell", cells[2])
                    blanket_filter = create_filter("cell", cells[3])
                    shield_filter = create_filter("cell", cells[4])
                    last_wall_filter = create_filter("cell", cells[5])
                elif q == 2:
                    first_wall_filter = create_filter("cell", cells[1])
                    multiplier_filter = create_filter("cell", cells[2])
                    blanket_filter = create_filter("cell", cells[3])
                    moderator_filter = create_filter("cell", cells[4])
                    last_wall_filter = create_filter("cell", cells[5])
                else:
                    first_wall_filter = create_filter("cell", cells[1])
                    multiplier_filter = create_filter("cell", cells[2])
                    blanket_filter = create_filter("cell", cells[3])
                    shield_filter = create_filter("cell", cells[5])
                    moderator_filter = create_filter("cell", cells[4])
                    last_wall_filter = create_filter("cell", cells[6])

                first_wall_tally = create_tally(
                    "first_wall_tally", [first_wall_filter], ["444"]
                )
                last_wall_tally = create_tally(
                    "last_wall_tally", [last_wall_filter], ["444"]
                )
                if q != 0:

                    multiplier_tally = create_tally(
                        "multiplier_tally", [multiplier_filter], ["444"]
                    )
                if q == 1 or q == 3:
                    shield_tally = create_tally(
                        "shield_tally", [shield_filter], ["444"]
                    )
                if q == 2 or q == 3:
                    moderator_tally = create_tally(
                        "moderator_tally", [moderator_filter], ["444"]
                    )
                blanket_tbr_tally = create_tally(
                    "blanket_tbr_tally", [blanket_filter], ["(n,Xt)"]
                )
                blanket_tbr_tally.nuclides = ["Li6", "Li7"]

                if q == 0:
                    tallies = openmc.Tallies(
                        [
                            first_wall_tally,
                            last_wall_tally,
                            blanket_tbr_tally,
                        ]
                    )
                elif q == 1:
                    tallies = openmc.Tallies(
                        [
                            first_wall_tally,
                            shield_tally,
                            multiplier_tally,
                            last_wall_tally,
                            blanket_tbr_tally,
                        ]
                    )
                elif q == 2:
                    tallies = openmc.Tallies(
                        [
                            first_wall_tally,
                            moderator_tally,
                            multiplier_tally,
                            last_wall_tally,
                            blanket_tbr_tally,
                        ]
                    )
                else:
                    tallies = openmc.Tallies(
                        [
                            first_wall_tally,
                            shield_tally,
                            moderator_tally,
                            multiplier_tally,
                            last_wall_tally,
                            blanket_tbr_tally,
                        ]
                    )

                model = openmc.model.Model(geometry, materials, settings, tallies)
                print(blanket_thickness)

                file_remover(".h5", dir_name)

                results_filename = model.run()
                results = openmc.StatePoint(results_filename)

                if q == 0:
                    first_wall_volume = sphere_vol(102, 100)
                    blanket_volume = sphere_vol(blanket_thickness + 102, 102)
                    last_wall_volume = sphere_vol(
                        blanket_thickness + 104, blanket_thickness + 102
                    )
                elif q == 1:
                    first_wall_volume = sphere_vol(102, 100)
                    multiplier_volume = sphere_vol(107, 102)
                    blanket_volume = sphere_vol(blanket_thickness + 107, 107)
                    shield_volume = sphere_vol(
                        blanket_thickness + 112, blanket_thickness + 107
                    )
                    last_wall_volume = sphere_vol(
                        blanket_thickness + 114, blanket_thickness + 112
                    )
                elif q == 2:
                    first_wall_volume = sphere_vol(102, 100)
                    multiplier_volume = sphere_vol(107, 102)
                    blanket_volume = sphere_vol(blanket_thickness + 107, 107)
                    moderator_volume = sphere_vol(
                        blanket_thickness + 112, blanket_thickness + 107
                    )
                    last_wall_volume = sphere_vol(
                        blanket_thickness + 114, blanket_thickness + 112
                    )
                else:
                    first_wall_volume = sphere_vol(102, 100)
                    multiplier_volume = sphere_vol(107, 102)
                    blanket_volume = sphere_vol(blanket_thickness + 107, 107)
                    moderator_volume = sphere_vol(
                        blanket_thickness + 112, blanket_thickness + 107
                    )
                    shield_volume = sphere_vol(
                        blanket_thickness + 117, blanket_thickness + 112
                    )
                    last_wall_volume = sphere_vol(
                        blanket_thickness + 119, blanket_thickness + 117
                    )

                first_wall_raw = results.get_tally(name="first_wall_tally")
                tbr_tally = results.get_tally(name="blanket_tbr_tally")
                last_wall_raw = results.get_tally(name="last_wall_tally")

                df_first_wall = first_wall_raw.get_pandas_dataframe()
                df_last_wall = last_wall_raw.get_pandas_dataframe()
                df_tbr = tbr_tally.get_pandas_dataframe()

                first_wall_damage_energy_in_ev = df_first_wall["mean"].sum()
                last_wall_damage_energy_in_ev = df_last_wall["mean"].sum()
                tbr_tally_result = df_tbr["mean"].sum()

                if q != 0:
                    multiplier_raw = results.get_tally(name="multiplier_tally")
                    df_multiplier = multiplier_raw.get_pandas_dataframe()
                    multiplier_damage_energy_in_ev = df_multiplier["mean"].sum()
                if q == 1 or q == 3:
                    shield_raw = results.get_tally(name="shield_tally")
                    df_shield = shield_raw.get_pandas_dataframe()
                    shield_damage_energy_in_ev = df_shield["mean"].sum()
                if q == 2 or q == 3:
                    moderator_raw = results.get_tally(name="moderator_tally")
                    df_moderator = moderator_raw.get_pandas_dataframe()
                    moderator_damage_energy_in_ev = df_moderator["mean"].sum()

                first_wall_dpa = calculate_dpa(
                    first_wall_damage_energy_in_ev, materials[0], first_wall_volume
                )

                if q == 0:
                    last_wall_dpa = calculate_dpa(
                        last_wall_damage_energy_in_ev, materials[2], last_wall_volume
                    )
                    temp_dpa.append(
                        [
                            first_wall_dpa,
                            last_wall_dpa,
                        ]
                    )
                elif q == 1:
                    multiplier_dpa = calculate_dpa(
                        multiplier_damage_energy_in_ev, materials[1], multiplier_volume
                    )
                    shield_dpa = calculate_dpa(
                        shield_damage_energy_in_ev, materials[3], shield_volume
                    )
                    last_wall_dpa = calculate_dpa(
                        last_wall_damage_energy_in_ev, materials[4], last_wall_volume
                    )
                    temp_dpa.append(
                        [first_wall_dpa, multiplier_dpa, shield_dpa, last_wall_dpa]
                    )
                elif q == 2:
                    multiplier_dpa = calculate_dpa(
                        multiplier_damage_energy_in_ev, materials[1], multiplier_volume
                    )
                    moderator_dpa = calculate_dpa(
                        moderator_damage_energy_in_ev, materials[3], moderator_volume
                    )
                    last_wall_dpa = calculate_dpa(
                        last_wall_damage_energy_in_ev, materials[4], last_wall_volume
                    )
                    temp_dpa.append(
                        [first_wall_dpa, multiplier_dpa, moderator_dpa, last_wall_dpa]
                    )
                elif q == 3:
                    multiplier_dpa = calculate_dpa(
                        multiplier_damage_energy_in_ev, materials[1], multiplier_volume
                    )
                    moderator_dpa = calculate_dpa(
                        moderator_damage_energy_in_ev, materials[3], moderator_volume
                    )
                    shield_dpa = calculate_dpa(
                        shield_damage_energy_in_ev, materials[4], shield_volume
                    )
                    last_wall_dpa = calculate_dpa(
                        last_wall_damage_energy_in_ev, materials[5], last_wall_volume
                    )
                    temp_dpa.append(
                        [
                            first_wall_dpa,
                            multiplier_dpa,
                            moderator_dpa,
                            shield_dpa,
                            last_wall_dpa,
                        ]
                    )
                with openmc.StatePoint("statepoint.10.h5") as sp:
                    leakage_mean = sp.global_tallies[3]["mean"]
                temp_leakage.append(leakage_mean)
                temp_tbr.append(tbr_tally_result)

            leakage_fraction.append(temp_leakage)
            dpa_values.append(temp_dpa)
            tbr_values.append(temp_tbr)

        blanket_sizes = list(range(10, 151, 5))
        x = 2
        if q == 0:
            x = 1
        materials = [material.name for material in mats[x]]

        colors = [
            "#FF0000",
            "#FFA500",
            "#008000",
            "#800080",
            "#0000FF",
            "#00FFFF",
            "#FF00FF",
            "#808080",
            "#008080",
            "#000080",
            "#000000",
            "#800000",
        ]

        first = ["first wall", "last wall"]
        second = ["first wall", "multiplier", "shield", "last wall"]
        third = ["first wall", "multiplier", "moderator", "last wall"]
        fourth = ["first wall", "multiplier", "moderator", "shield", "last wall"]
        layers = [first, second, third, fourth]

        print(dpa_values)

        for a in range(len(layers[q])):
            # Extract DPA values for the current layer
            layer_name = layers[q][a]
            dpa_layer = []
            for b in range(len(dpa_values)):
                temp = []
                for c in range(len(dpa_values[b])):
                    temp.append(dpa_values[b][c][a])
                dpa_layer.append(temp)

            dpa_layer_transposed = list(map(list, zip(*dpa_layer)))
            df_layer = pd.DataFrame(
                dpa_layer_transposed, columns=blanket_sizes, index=materials
            )

            # Create a scientific DPA graph using Plotly
            fig_dpa = go.Figure()

            for idx, material in enumerate(df_layer.index):
                color = colors[idx % len(colors)]
                fig_dpa.add_trace(
                    go.Scatter(
                        x=df_layer.columns,
                        y=df_layer.loc[material],
                        mode="lines",
                        name=f"{material}",
                        line=dict(color=color),
                    )
                )

            # Set the x and y axis to linear scale
            fig_dpa.update_xaxes(
                type="linear", title="<b>THICKNESS OF BLANKET (CM)</b>"
            )
            fig_dpa.update_yaxes(
                type="linear", title="<b>DISPLACEMENT PER ATOM (DPA)</b>"
            )

            # Remove the title
            fig_dpa.update_layout(
                legend=dict(
                    font=dict(size=10),
                    borderwidth=1,
                    bordercolor="black",
                    title=dict(text="<b>LEGEND</b>", font=dict(size=12)),
                ),  # Add legend boundary and title
                showlegend=True,
                xaxis_showgrid=True,
                yaxis_showgrid=True,
                margin=dict(l=20, r=20, t=50, b=20),  # Add slight margin
                title=f"<b>DISPLACEMENT PER ATOM FOR {layer_name.upper()}</b>",
                title_x=0.5,
            )

            # Save the figure
            pio.write_image(fig_dpa, f"graphs/{q}/DPA_{layer_name}.png")

            # Clear the figure
            fig_dpa = None

        tbr_values_transposed = list(map(list, zip(*tbr_values)))
        leakage_values_transposed = list(map(list, zip(*leakage_fraction)))

        # Create a DataFrame from tbr_values and leakage_values
        df_tbr = pd.DataFrame(
            tbr_values_transposed, columns=blanket_sizes, index=materials
        )
        df_leakage = pd.DataFrame(
            leakage_values_transposed, columns=blanket_sizes, index=materials
        )

        fig_tbr = go.Figure()

        colors_used = []
        for idx, material in enumerate(df_tbr.index):
            color = colors[idx % len(colors)]  # Use a different color for each material
            colors_used.append(color)
            fig_tbr.add_trace(
                go.Scatter(
                    x=df_tbr.columns,
                    y=df_tbr.loc[material],
                    mode="lines",
                    name=f"{material} TBR",
                    line=dict(color=color),
                )
            )

        fig_tbr.update_xaxes(type="linear")
        fig_tbr.update_yaxes(type="linear")

        # Set the x and y axis labels
        fig_tbr.update_layout(
            xaxis_title="<b>BLANKET THICKNESS (CM)</b>",
            yaxis_title="<b>TBR</b>",
            title="<b>TBR VS. BLANKET THICKNESS</b>",
            title_x=0.5,  # Center the title
            margin=dict(l=20, r=20, b=20, t=50),  # Add slight margin
            legend=dict(
                font=dict(size=10),
                borderwidth=1,
                bordercolor="black",
                title=dict(text="<b>LEGEND</b>", font=dict(size=12)),
            ),  # Add legend boundary and title
        )
        pio.write_image(fig_tbr, f"graphs/{q}/TBR.png")

        for idx, material in enumerate(df_tbr.index):
            fig_tbr.add_trace(
                go.Scatter(
                    x=df_leakage.columns,
                    y=df_leakage.loc[material],
                    mode="lines",
                    name=f"{material} Leakage",
                    line=dict(dash="dash", color=colors_used[idx]),
                )
            )

        # Set the x and y axis to linear scale
        fig_tbr.update_xaxes(type="linear")
        fig_tbr.update_yaxes(type="linear")

        # Set the x and y axis labels
        fig_tbr.update_layout(
            xaxis_title="<b>BLANKET THICKNESS (CM)</b>",
            yaxis_title="<b>TBR / LEAKAGE</b>",
            title="<b>TBR AND LEAKAGE VS. BLANKET THICKNESS</b>",
            title_x=0.5,  # Center the title
            margin=dict(l=20, r=20, b=20, t=50),  # Add slight margin
            legend=dict(
                font=dict(size=10),
                borderwidth=1,
                bordercolor="black",
                title=dict(text="<b>LEGEND</b>", font=dict(size=12)),
            ),  # Add legend boundary and title
        )

        # Save the modified PNG file
        pio.write_image(fig_tbr, f"graphs/{q}/TBR_and_Leakage.png")

        fig_tbr = go.Figure()

        for idx, material in enumerate(df_tbr.index):
            fig_tbr.add_trace(
                go.Scatter(
                    x=df_leakage.columns,
                    y=df_leakage.loc[material],
                    mode="lines",
                    name=f"{material} Leakage",
                    line=dict(color=colors_used[idx]),
                )
            )

        # Set the x and y axis to linear scale
        fig_tbr.update_xaxes(type="linear")
        fig_tbr.update_yaxes(type="linear")

        # Set the x and y axis labels
        fig_tbr.update_layout(
            xaxis_title="<b>BLANKET THICKNESS (CM)</b>",
            yaxis_title="<b>LEAKAGE</b>",
            title="<b>LEAKAGE VS. BLANKET THICKNESS</b>",
            title_x=0.5,  # Center the title
            margin=dict(l=20, r=20, b=20, t=50),  # Add slight margin
            legend=dict(
                font=dict(size=10),
                borderwidth=1,
                bordercolor="black",
                title=dict(text="<b>LEGEND</b>", font=dict(size=12)),
            ),  # Add legend boundary and title
        )

        # Save the modified PNG file
        pio.write_image(fig_tbr, f"graphs/{q}/Leakage.png")


if __name__ == "__main__":
    dpa_and_tbr()
