"""
Main script to plot flow graphs from ROOT files using JPyPlotRatio.
Reads v2_2_eg1.0 and v3_2_eg1.0 vs centrality for OO and NeNe systems, comparing AC vs WS.
"""

import sys
sys.path.append("JPyPlotRatio")
import numpy as np
import ROOT
import JPyPlotRatio
import array

# Plot parameters for theory curves (light colors with transparency)
# Order matches fileIndices: [WS, AC, ACcompact, ACloose]
# AC variants use blue color scheme (violet, blue, teal), WS uses different color
plotParameters = [
    {'plotType':'theory','color':'#FFA07A','linecolor':'#FFA07A', 'linestyle':'--','linewidth':2,'alpha':0.6,'labelLegendId':0},    # WS - light salmon (different color)
    {'plotType':'theory','color':'#DDA0DD','linecolor':'#DDA0DD','linestyle':'-','linewidth':2,'alpha':0.6,'labelLegendId':0},    # AC - plum (light violet)
    {'plotType':'theory','color':'#4682B4','linecolor':'#4682B4','linestyle':':','linewidth':2,'alpha':0.6,'labelLegendId':0},    # ACcompact - steel blue (darker blue)
    {'plotType':'theory','color':'#AFEEEE','linecolor':'#AFEEEE','linestyle':'-.','linewidth':2,'alpha':0.6,'labelLegendId':0},    # ACloose - pale turquoise (light teal)
]

# File paths: [OO WS, OO AC, OO ACcompact, OO ACloose, NeNe WS, NeNe AC, NeNe ACcompact, NeNe ACloose]
filenames = [
    "flow/vncn_OOWS.root",           # Index 0: OO WS
    "flow/vncn_OOAC.root",            # Index 1: OO AC
    "flow/vncn_OOACcompact.root",     # Index 2: OO ACcompact
    "flow/vncn_OOACloose.root",       # Index 3: OO ACloose
    "flow/vncn_NeNeWS.root",          # Index 4: NeNe WS
    "flow/vncn_NeNeAC.root",          # Index 5: NeNe AC
    "flow/vncn_NeNeACcompact.root",   # Index 6: NeNe ACcompact
    "flow/vncn_NeNeACloose.root"      # Index 7: NeNe ACloose
]

# Graph names
graphNames = {
    'v2': 'v2_2_eg1.0',
    'v3': 'v3_2_eg1.0'
}

# Centrality bins (midpoints of 5% bins)
centrality_bins = np.array([2.5, 7.5, 12.5, 17.5, 22.5, 27.5, 32.5, 37.5, 42.5, 47.5, 52.5, 57.5])

# ALICE data 2509.06428
#errors are missing
v2_OO_data_alice = np.array([0.060, 0.062, 0.065, 0.068, 0.070, 0.071, 0.070, 0.069, 0.068, 0.067, 0.066, 0.065])
v2_NeNe_data_alice = np.array([0.060, 0.061, 0.064, 0.067, 0.069, 0.070, 0.069, 0.068, 0.067, 0.066, 0.065, 0.064])
v3_OO_data_alice = np.array([0.023, 0.022, 0.020, 0.018, 0.016, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])
v3_NeNe_data_alice = np.array([0.022, 0.021, 0.019, 0.017, 0.015, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010, 0.009])

# CMS data 2510.02580
#errors are missing
v2_OO_data_cms = np.array([0.060, 0.062, 0.065, 0.070, 0.078, 0.080, 0.079, 0.077, 0.075, 0.073, 0.071, 0.070])
v2_NeNe_data_cms = np.array([0.060, 0.062, 0.065, 0.070, 0.078, 0.080, 0.079, 0.077, 0.075, 0.073, 0.071, 0.070])
v3_OO_data_cms = np.array([0.025, 0.024, 0.022, 0.020, 0.018, 0.016, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])
v3_NeNe_data_cms = np.array([0.025, 0.024, 0.022, 0.020, 0.018, 0.016, 0.014, 0.013, 0.012, 0.011, 0.010, 0.010])

# Store data in dictionaries for easy access
alice_data = {
    ('v2', 'OO'): (centrality_bins, v2_OO_data_alice),
    ('v2', 'NeNe'): (centrality_bins, v2_NeNe_data_alice),
    ('v3', 'OO'): (centrality_bins, v3_OO_data_alice),
    ('v3', 'NeNe'): (centrality_bins, v3_NeNe_data_alice)
}

cms_data = {
    ('v2', 'OO'): (centrality_bins, v2_OO_data_cms),
    ('v2', 'NeNe'): (centrality_bins, v2_NeNe_data_cms),
    ('v3', 'OO'): (centrality_bins, v3_OO_data_cms),
    ('v3', 'NeNe'): (centrality_bins, v3_NeNe_data_cms)
}

# Panel layout: 2 rows (v2, v3) x 2 columns (OO, NeNe) = 4 panels
# Panel indices: 0=v2_OO, 1=v2_NeNe, 2=v3_OO, 3=v3_NeNe
plot = JPyPlotRatio.JPyPlotRatio(
    panels=(2, 2), 
    panelsize=(10, 8), 
    logScale=0,
    axisLabelSize=20, 
    tickLabelSize=16, 
    rowBounds={0: (0.02, 0.12), 1: (0.0, 0.05)},
    colBounds={0: (-5.0, 90.0), 1: (-5.0, 90.0)},
    xlabel="Centrality (%)",
    ylabel="$v_n$", 
    legendPanel={0: 0, 1: 1},
    disableRatio={},
    legendLoc={0: (0.6, 0.75), 1: (0.6, 0.75)}, 
    ratioBounds={0: (0.7, 1.3), 1: (0.2, 1.8)},
    legendSize=12
)

legends = ["WS", "AC (Default)", "AC (Compact)", "AC (Loose)"]

# Variables and systems to plot
variables = ['v2', 'v3']
systems = ['OO', 'NeNe']

legendIndex = 0
for variable in variables:
    for system in systems:
        graphName = graphNames[variable]
        
        print(f"Plotting {variable} for {system} (panel {legendIndex})...")
        
        # Determine file indices matching filenames array:
        # OO: indices [0, 1, 2, 3] = [WS, AC, ACcompact, ACloose]
        # NeNe: indices [4, 5, 6, 7] = [WS, AC, ACcompact, ACloose]
        if system == 'OO':
            fileIndices = [0, 1, 2, 3]  # WS, AC, ACcompact, ACloose
        else:  # NeNe
            fileIndices = [4, 5, 6, 7]  # WS, AC, ACcompact, ACloose
        
        p = []
        for colorIndex, fileIdx in enumerate(fileIndices):
            graph = ROOT.TFile(filenames[fileIdx], "read").Get(graphName)
            if not graph:
                print(f"Error: Graph {graphName} not found in {filenames[fileIdx]}")
                continue
            p.append(plot.Add(legendIndex, graph, label=legends[colorIndex], **plotParameters[colorIndex]))
        
        # Store references to theory graphs (WS=p[0], AC=p[1], ACcompact=p[2], ACloose=p[3])
        theory_ws = p[0] if len(p) > 0 else None
        theory_ac = p[1] if len(p) > 1 else None
        
        # Add ALICE data points as TGraphErrors (data plotType requires error bars)
        alice_graph_ref = None
        cms_graph_ref = None
        
        if (variable, system) in alice_data:
            x_data, y_data = alice_data[(variable, system)]
            # Create TGraphErrors with zero errors (data plotType requires error bars)
            x_arr = array.array('d', x_data.tolist())
            y_arr = array.array('d', y_data.tolist())
            ex_arr = array.array('d', [0.0] * len(x_data))  # Zero x errors
            ey_arr = array.array('d', [0.0] * len(y_data))  # Zero y errors
            alice_graph = ROOT.TGraphErrors(len(x_data), x_arr, y_arr, ex_arr, ey_arr)
            
            # Verify graph was created
            if alice_graph.GetN() == 0:
                print(f"Warning: ALICE graph for {variable} {system} has no points")
                continue
            
            # Plot parameters for ALICE data (v2=red, v3=blue)
            if variable == 'v2':
                alice_color = '#d62728'  # Red for v2
                alice_edge_color = '#b71c1c'  # Darker red edge
            else:  # v3
                alice_color = '#1f77b4'  # Blue for v3
                alice_edge_color = '#1565c0'  # Darker blue edge
            
            alice_params = {
                'plotType': 'data',
                'color': alice_color,
                'marker': 'o',  # Circle for both OO and NeNe
                'markersize': 10,
                'markeredgecolor': alice_edge_color,
                'markeredgewidth': 2.0,
                'markerfacecolor': alice_color,
                'alpha': 0.6,
                'labelLegendId': 1
            }
            
            # Add to plot (this will add to legend automatically)
            alice_graph_ref = plot.Add(legendIndex, alice_graph, label=f"ALICE ${variable.replace('v', 'v_')}$", **alice_params)
        
        # Add CMS data points as TGraphErrors
        if (variable, system) in cms_data:
            x_data, y_data = cms_data[(variable, system)]
            # Create TGraphErrors with zero errors
            x_arr = array.array('d', x_data.tolist())
            y_arr = array.array('d', y_data.tolist())
            ex_arr = array.array('d', [0.0] * len(x_data))
            ey_arr = array.array('d', [0.0] * len(y_data))
            cms_graph = ROOT.TGraphErrors(len(x_data), x_arr, y_arr, ex_arr, ey_arr)
            
            if cms_graph.GetN() == 0:
                print(f"Warning: CMS graph for {variable} {system} has no points")
            else:
                # Plot parameters for CMS data (same colors as ALICE but different markers)
                if variable == 'v2':
                    cms_color = '#d62728'  # Red for v2 (same as ALICE)
                    cms_edge_color = '#b71c1c'  # Darker red edge (same as ALICE)
                    cms_marker = 'P'  # Big plus for v2
                else:  # v3
                    cms_color = '#1f77b4'  # Blue for v3 (same as ALICE)
                    cms_edge_color = '#1565c0'  # Darker blue edge (same as ALICE)
                    cms_marker = 'X'  # Big X for v3
                
                cms_params = {
                    'plotType': 'data',
                    'color': cms_color,
                    'marker': cms_marker,
                    'markersize': 12,  # Big markers
                    'markeredgecolor': cms_edge_color,
                    'markeredgewidth': 2.5,
                    'markerfacecolor': 'none',  # Open markers (no fill)
                    'alpha': 0.7,
                    'labelLegendId': 1
                }
                
                # Add to plot
                cms_graph_ref = plot.Add(legendIndex, cms_graph, label=f"CMS ${variable.replace('v', 'v_')}$", **cms_params)
        
        # Add ratios: data over MC (only WS ratios)
        if theory_ws is not None:
            if alice_graph_ref is not None:
                plot.Ratio(alice_graph_ref, theory_ws)  # ALICE/WS
            if cms_graph_ref is not None:
                plot.Ratio(cms_graph_ref, theory_ws)  # CMS/WS
        
        # Add text labels (bigger text for all panels)
        ax = plot.GetAxes(legendIndex)
        ax.text(0.05, 0.9, system, size=18, transform=ax.transAxes)
        
        # Add kinematics information (only on first panel)
        if legendIndex == 0:  # v2_OO panel only
            ax.text(0.05, 0.15, "ALICE: $0.2 < p_T < 3.0$ GeV/$c$, $|\\eta| < 0.8$", 
                   size=12, transform=ax.transAxes, color='#666666')
            ax.text(0.05, 0.10, "CMS: $0.3 < p_T < 3.0$ GeV/$c$, $|\\eta| < 2.4$, $|\\Delta\\eta| > 2$", 
                   size=12, transform=ax.transAxes, color='#666666')
        
        legendIndex += 1

# Add row labels (bigger text for all)
plot.GetAxes(0).text(0.05, 0.8, "$v_2$", size=20, transform=plot.GetAxes(0).transAxes, weight='bold')
plot.GetAxes(1).text(0.05, 0.8, "$v_2$", size=20, transform=plot.GetAxes(1).transAxes, weight='bold')
plot.GetAxes(2).text(0.05, 0.8, "$v_3$", size=20, transform=plot.GetAxes(2).transAxes, weight='bold')
plot.GetAxes(3).text(0.05, 0.8, "$v_3$", size=20, transform=plot.GetAxes(3).transAxes, weight='bold')

plot.Plot()
plot.Save('flow_comparison.pdf')
print("Saved: flow_comparison.pdf")
