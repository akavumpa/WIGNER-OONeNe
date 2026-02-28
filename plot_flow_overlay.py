"""
Overlay plot: ALICE, ATLAS, and CMS flow data only. No theory.
Ratio = data / ALICE (so ATLAS/ALICE and CMS/ALICE).
Axis range as ATLAS. Different colors and markers per experiment.
Output: flow_comparison_overlay.pdf
"""

import sys
sys.path.append("JPyPlotRatio")
import numpy as np
import ROOT
import JPyPlotRatio
import array

centrality_bins = np.array([2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5])

# ALICE
v2_OO_data_alice = np.array([0.060, 0.062, 0.065, 0.068, 0.070, 0.071, 0.070, 0.069, 0.068, 0.067, 0.066, 0.065])
v2_NeNe_data_alice = np.array([0.060, 0.061, 0.064, 0.067, 0.069, 0.070, 0.069, 0.068, 0.067, 0.066, 0.065, 0.064])
v3_OO_data_alice = np.array([0.023, 0.022, 0.020, 0.018, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])
v3_NeNe_data_alice = np.array([0.022, 0.021, 0.019, 0.017, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010, 0.009])

# CMS
v2_OO_data_cms = np.array([0.060, 0.062, 0.065, 0.070, 0.078, 0.080, 0.079, 0.077, 0.075, 0.073, 0.071, 0.070])
v2_NeNe_data_cms = np.array([0.060, 0.062, 0.065, 0.070, 0.078, 0.080, 0.079, 0.077, 0.075, 0.073, 0.071, 0.070])
v3_OO_data_cms = np.array([0.025, 0.024, 0.022, 0.020, 0.018, 0.016, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])
v3_NeNe_data_cms = np.array([0.025, 0.024, 0.022, 0.020, 0.018, 0.016, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])

# ATLAS
v2_OO_data_atlas = np.array([0.082, 0.088, 0.092, 0.094, 0.095, 0.094, 0.092, 0.090, 0.082, 0.075, 0.070, 0.068])
v2_NeNe_data_atlas = np.array([0.085, 0.090, 0.093, 0.095, 0.096, 0.095, 0.093, 0.091, 0.084, 0.076, 0.071, 0.069])
v3_OO_data_atlas = np.array([0.038, 0.036, 0.034, 0.032, 0.030, 0.029, 0.027, 0.026, 0.024, 0.022, 0.020, 0.019])
v3_NeNe_data_atlas = np.array([0.040, 0.038, 0.036, 0.034, 0.032, 0.031, 0.029, 0.028, 0.026, 0.024, 0.022, 0.021])

alice_data = {('v2', 'OO'): (centrality_bins, v2_OO_data_alice), ('v2', 'NeNe'): (centrality_bins, v2_NeNe_data_alice),
              ('v3', 'OO'): (centrality_bins, v3_OO_data_alice), ('v3', 'NeNe'): (centrality_bins, v3_NeNe_data_alice)}
cms_data = {('v2', 'OO'): (centrality_bins, v2_OO_data_cms), ('v2', 'NeNe'): (centrality_bins, v2_NeNe_data_cms),
            ('v3', 'OO'): (centrality_bins, v3_OO_data_cms), ('v3', 'NeNe'): (centrality_bins, v3_NeNe_data_cms)}
atlas_data = {('v2', 'OO'): (centrality_bins, v2_OO_data_atlas), ('v2', 'NeNe'): (centrality_bins, v2_NeNe_data_atlas),
              ('v3', 'OO'): (centrality_bins, v3_OO_data_atlas), ('v3', 'NeNe'): (centrality_bins, v3_NeNe_data_atlas)}

# Data style: (exp_name, data_dict, color, edge_color, marker_v2, marker_v3, filled)
# ALICE first (ratio reference), then ATLAS, CMS
experiments = [
    ('ALICE', alice_data, '#d62728', '#b71c1c', 'o', 's', True),
    ('ATLAS', atlas_data, '#1976D2', '#0D47A1', '^', 'v', True),
    ('CMS', cms_data, '#388E3C', '#1B5E20', 'P', 'X', False),
]

variables = ['v2', 'v3']
systems = ['OO', 'NeNe']

# ATLAS axis range
row_bounds = {0: (0.035, 0.12), 1: (0.0, 0.06)}

plot = JPyPlotRatio.JPyPlotRatio(
    panels=(2, 2),
    panelsize=(10, 8),
    logScale=0,
    axisLabelSize=20,
    tickLabelSize=16,
    rowBounds=row_bounds,
    colBounds={0: (-5.0, 90.0), 1: (-5.0, 90.0)},
    xlabel="Centrality (%)",
    ylabel="$v_n$",
    legendPanel={0: 0, 1: 0, 2:1},
    disableRatio={},
    legendLoc={0: (0.4, 0.85), 1: (0.7, 0.85), 2: (0.7, 0.85)},
    ratioBounds={0: (0.5, 1.8), 1: (0.6, 2.6)},
    legendSize=16,
)

legendIndex = 0
for variable in variables:
    for system in systems:
        alice_ref = None
        for exp_name, data_dict, color, edge_color, m2, m3, filled in experiments:
            if (variable, system) not in data_dict:
                continue
            x_data, y_data = data_dict[(variable, system)]
            x_arr = array.array('d', x_data.tolist())
            y_arr = array.array('d', y_data.tolist())
            ex_arr = array.array('d', [0.0] * len(x_data))
            ey_arr = array.array('d', [0.0] * len(y_data))
            g = ROOT.TGraphErrors(len(x_data), x_arr, y_arr, ex_arr, ey_arr)
            if g.GetN() == 0:
                continue
            marker = m2 if variable == 'v2' else m3
            params = {
                'plotType': 'data', 'color': color, 'marker': marker,
                'markersize': 10, 'markeredgecolor': edge_color, 'markeredgewidth': 2.0,
                'markerfacecolor': color if filled else 'none',
                'alpha': 0.7, 'labelLegendId': 0 if exp_name == 'ALICE' else 1 if exp_name == 'ATLAS' else 2,
            }
            if not filled:
                params['markersize'] = 12
                params['markeredgewidth'] = 2.5
            ref = plot.Add(legendIndex, g, label=f"{exp_name} ${variable.replace('v', 'v_')}$", **params)
            if exp_name == 'ALICE':
                alice_ref = ref
            else:
                # Ratio = data / ALICE
                if alice_ref is not None:
                    plot.Ratio(ref, alice_ref)

        ax = plot.GetAxes(legendIndex)
        ax.text(0.05, 0.9, system, size=18, transform=ax.transAxes)
        if legendIndex == 0:
            ax.text(0.05, 0.12, "ALICE, ATLAS, CMS | ratio to ALICE", size=11, transform=ax.transAxes, color='#666666')
        legendIndex += 1

for i in range(4):
    plot.GetAxes(i).text(0.05, 0.8, "$v_2$" if i < 2 else "$v_3$", size=20, transform=plot.GetAxes(i).transAxes, weight='bold')

plot.Plot()
plot.Save('flow_comparison_overlay.pdf')
print("Saved: flow_comparison_overlay.pdf")
