"""
Plot flow v2 and v3 vs centrality for OO and NeNe: ALICE data vs ALICE-like theory,
and CMS data vs CMS-like theory. Produces two separate plots.
Theory: flow/vn_Centrality/ALICE/ and flow/vn_Centrality/CMS/ (vn_OO.root, vn_NeNe.root,
vn_OOVary.root, vn_NeNeVary.root) with graph names v2_OOAC, v2_OOWS, v2_OOCompact, v2_OOLoose, etc.
"""

import sys
sys.path.append("JPyPlotRatio")
import numpy as np
import ROOT
import JPyPlotRatio
import array

# Plot parameters: WS = distinct color; AC family = blue, green (dark but not so), light purple
plotParameters = [
    {'plotType': 'theory', 'color': '#E07A5F', 'linecolor': '#E07A5F', 'linestyle': '--', 'linewidth': 2, 'alpha': 0.7, 'labelLegendId': 0},   # WS (terracotta - distinct)
    {'plotType': 'theory', 'color': '#2196F3', 'linecolor': '#2196F3', 'linestyle': '-', 'linewidth': 2, 'alpha': 0.6, 'labelLegendId': 0},    # AC (blue)
    {'plotType': 'theory', 'color': '#388E3C', 'linecolor': '#388E3C', 'linestyle': ':', 'linewidth': 2, 'alpha': 0.6, 'labelLegendId': 0},    # AC Compact (green, dark but not so)
    {'plotType': 'theory', 'color': '#B39DDB', 'linecolor': '#B39DDB', 'linestyle': '-.', 'linewidth': 2, 'alpha': 0.6, 'labelLegendId': 0},   # AC Loose (light purple)
]

legends_theory = ["WS", "AC (Default)", "AC (Compact)", "AC (Loose)"]

# Graph names: vn_OO.root/vn_NeNe.root use v2_OOAC, v2_OOWS, v3_OOAC, v3_OOWS; v2_NeNeAC, v2_NeNeWS, etc.
# vn_OOVary.root/vn_NeNeVary.root use v2_OOCompact, v2_OOLoose, v3_OOCompact, v3_OOLoose; v2_NeNeCompact, etc.
def get_graph_names(variable, system):
    return {
        'WS': f"{variable}_{system}WS",      # v2_OOWS, v2_NeNeWS, ...
        'AC': f"{variable}_{system}AC",      # v2_OOAC, v2_NeNeAC, ...
        'Compact': f"{variable}_{system}Compact",
        'Loose': f"{variable}_{system}Loose",
    }

# Centrality bins (midpoints of 5% bins)
centrality_bins = np.array([2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5])

# ALICE data 2509.06428
v2_OO_data_alice = np.array([0.060, 0.062, 0.065, 0.068, 0.070, 0.071, 0.070, 0.069, 0.068, 0.067, 0.066, 0.065])
v2_NeNe_data_alice = np.array([0.060, 0.061, 0.064, 0.067, 0.069, 0.070, 0.069, 0.068, 0.067, 0.066, 0.065, 0.064])
v3_OO_data_alice = np.array([0.023, 0.022, 0.020, 0.018, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])
v3_NeNe_data_alice = np.array([0.022, 0.021, 0.019, 0.017, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010, 0.009])

# CMS data 2510.02580
v2_OO_data_cms = np.array([0.060, 0.062, 0.065, 0.070, 0.078, 0.080, 0.079, 0.077, 0.075, 0.073, 0.071, 0.070])
v2_NeNe_data_cms = np.array([0.060, 0.062, 0.065, 0.070, 0.078, 0.080, 0.079, 0.077, 0.075, 0.073, 0.071, 0.070])
v3_OO_data_cms = np.array([0.025, 0.024, 0.022, 0.020, 0.018, 0.016, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])
v3_NeNe_data_cms = np.array([0.025, 0.024, 0.022, 0.020, 0.018, 0.016, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])

alice_data = {
    ('v2', 'OO'): (centrality_bins, v2_OO_data_alice),
    ('v2', 'NeNe'): (centrality_bins, v2_NeNe_data_alice),
    ('v3', 'OO'): (centrality_bins, v3_OO_data_alice),
    ('v3', 'NeNe'): (centrality_bins, v3_NeNe_data_alice),
}

cms_data = {
    ('v2', 'OO'): (centrality_bins, v2_OO_data_cms),
    ('v2', 'NeNe'): (centrality_bins, v2_NeNe_data_cms),
    ('v3', 'OO'): (centrality_bins, v3_OO_data_cms),
    ('v3', 'NeNe'): (centrality_bins, v3_NeNe_data_cms),
}

variables = ['v2', 'v3']
systems = ['OO', 'NeNe']


def make_flow_plot(exp_name, data_dict, theory_dir, output_pdf, use_cms_markers=False):
    """
    Build one 2x2 flow plot: data_dict vs theory from theory_dir.
    exp_name: 'ALICE' or 'CMS'
    use_cms_markers: if True, use P for v2 and X for v3 data; else circles.
    """
    plot = JPyPlotRatio.JPyPlotRatio(
        panels=(2, 2),
        panelsize=(10, 8),
        logScale=0,
        axisLabelSize=20,
        tickLabelSize=16,
        rowBounds={0: (0.035, 0.12), 1: (0.0, 0.05)},
        colBounds={0: (-5.0, 90.0), 1: (-5.0, 90.0)},
        xlabel="Centrality (%)",
        ylabel="$v_n$",
        legendPanel={0: 0, 1: 1},
        disableRatio={},
        legendLoc={0: (0.6, 0.8), 1: (0.6, 0.8)},
        ratioBounds={0: (0.5, 1.5), 1: (0.6, 2.4)},
        legendSize=16,
    )

    legendIndex = 0
    for variable in variables:
        for system in systems:
            names = get_graph_names(variable, system)
            # Files: vn_OO.root, vn_NeNe.root for WS, AC; vn_OOVary.root, vn_NeNeVary.root for Compact, Loose
            file_main = f"{theory_dir}/vn_{system}.root"
            file_vary = f"{theory_dir}/vn_{system}Vary.root"

            theory_refs = []
            for i, (key, graph_name) in enumerate([('WS', names['WS']), ('AC', names['AC'])]):
                f = ROOT.TFile.Open(file_main, "read")
                if not f:
                    print(f"Error: Could not open {file_main}")
                    continue
                g = f.Get(graph_name)
                if not g:
                    print(f"Error: Graph {graph_name} not found in {file_main}")
                    f.Close()
                    continue
                ref = plot.Add(legendIndex, g, label=legends_theory[i], **plotParameters[i])
                theory_refs.append(ref)
                f.Close()

            for i, (key, graph_name) in enumerate([('Compact', names['Compact']), ('Loose', names['Loose'])]):
                f = ROOT.TFile.Open(file_vary, "read")
                if not f:
                    print(f"Error: Could not open {file_vary}")
                    continue
                g = f.Get(graph_name)
                if not g:
                    print(f"Error: Graph {graph_name} not found in {file_vary}")
                    f.Close()
                    continue
                ref = plot.Add(legendIndex, g, label=legends_theory[2 + i], **plotParameters[2 + i])
                theory_refs.append(ref)
                f.Close()

            # Data points
            data_ref = None
            if (variable, system) in data_dict:
                x_data, y_data = data_dict[(variable, system)]
                x_arr = array.array('d', x_data.tolist())
                y_arr = array.array('d', y_data.tolist())
                ex_arr = array.array('d', [0.0] * len(x_data))
                ey_arr = array.array('d', [0.0] * len(y_data))
                data_graph = ROOT.TGraphErrors(len(x_data), x_arr, y_arr, ex_arr, ey_arr)
                if data_graph.GetN() == 0:
                    print(f"Warning: Data graph for {variable} {system} has no points")
                else:
                    if variable == 'v2':
                        color = '#d62728'
                        edge_color = '#b71c1c'
                    else:
                        color = '#1f77b4'
                        edge_color = '#1565c0'
                    if use_cms_markers:
                        marker = 'P' if variable == 'v2' else 'X'
                        data_params = {
                            'plotType': 'data', 'color': color, 'marker': marker,
                            'markersize': 12, 'markeredgecolor': edge_color, 'markeredgewidth': 2.5,
                            'markerfacecolor': 'none', 'alpha': 0.7, 'labelLegendId': 1,
                        }
                    else:
                        # ALICE: v2 = circle (sphere), v3 = square
                        alice_marker = 'o' if variable == 'v2' else 's'
                        data_params = {
                            'plotType': 'data', 'color': color, 'marker': alice_marker,
                            'markersize': 10, 'markeredgecolor': edge_color, 'markeredgewidth': 2.0,
                            'markerfacecolor': color, 'alpha': 0.6, 'labelLegendId': 1,
                        }
                    data_ref = plot.Add(legendIndex, data_graph, label=f"{exp_name} ${variable.replace('v', 'v_')}$", **data_params)
                    # Ratio = model / data for all models
                    if data_ref is not None:
                        for theory_ref in theory_refs:
                            plot.Ratio(theory_ref, data_ref)

            ax = plot.GetAxes(legendIndex)
            ax.text(0.05, 0.9, system, size=18, transform=ax.transAxes)
            if legendIndex == 0:
                if exp_name == "ALICE":
                    ax.text(0.05, 0.15, f"{exp_name}: $0.2 < p_T < 3.0$ GeV/$c$, $|\\eta| < 0.8$", size=12, transform=ax.transAxes, color='#666666')
                else:
                    ax.text(0.05, 0.15, f"{exp_name}: $0.3 < p_T < 3.0$ GeV/$c$, $|\\eta| < 2.4$, $|\\Delta\\eta| > 2$", size=12, transform=ax.transAxes, color='#666666')
            legendIndex += 1

    for i in range(4):
        plot.GetAxes(i).text(0.05, 0.8, "$v_2$" if i < 2 else "$v_3$", size=20, transform=plot.GetAxes(i).transAxes, weight='bold')

    plot.Plot()
    plot.Save(output_pdf)
    print(f"Saved: {output_pdf}")


# ALICE: ALICE data vs ALICE-like theory
make_flow_plot(
    exp_name="ALICE",
    data_dict=alice_data,
    theory_dir="flow/vn_Centrality/ALICE",
    output_pdf="flow_comparison_ALICE.pdf",
    use_cms_markers=False,
)

# CMS: CMS data vs CMS-like theory
make_flow_plot(
    exp_name="CMS",
    data_dict=cms_data,
    theory_dir="flow/vn_Centrality/CMS",
    output_pdf="flow_comparison_CMS.pdf",
    use_cms_markers=True,
)
