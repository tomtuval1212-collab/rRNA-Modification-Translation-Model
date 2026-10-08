import json
import os
import random
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.lines import Line2D


# Set working paths
script_dir = os.path.dirname(os.path.abspath(__file__)) if '__file__' in globals() else os.getcwd() #gets the folder path where your main.py file is
data_dir = os.path.join(script_dir) #sets the input data directory path to the main project folder
output_dir = os.path.join(script_dir, 'output') #sets the output file (plots) to be in a folder named output
os.makedirs(output_dir, exist_ok=True) #creates the output folder unless it already exists



# -------------------------------------------------------------------
# 1. LOAD DATA FILES & BUILD REGISTRIES
# -------------------------------------------------------------------
#creating path strings do data files
viral_json_path = os.path.join(data_dir, 'viral_frameshift_sites.json')
viral_fasta_path = os.path.join(data_dir, 'viral_sequences.fasta')
params_path = os.path.join(data_dir, 'ribosome_params.json')
rmbase_path = os.path.join(data_dir, 'rrna_modifications_curated.csv')

#Loading JSON configuration files as dictionaries and CSV data into a DataFrame
try:
    with open(viral_json_path, 'r', encoding='utf-8') as f:
        viral_metadata = json.load(f)

    with open(params_path, 'r', encoding='utf-8') as f:
        ribo_params = json.load(f)

    df_rrna_mods = pd.read_csv(rmbase_path)

    df_dwell = pd.read_csv(os.path.join(data_dir, 'riboseq_dwell_times.csv')) #loades codons and their dwell times into a dataframe
except FileNotFoundError as e:
    raise SystemExit(f"Error loading required project files: {e}")

dwell_map = dict(zip(df_dwell['codon'], df_dwell['wt_dwell_time_sec'])) #takes codon and its dwell time from dataframe and creates a dictionary for a faster lookup
FOOTPRINT_NTS = int(ribo_params.get('ribosome_footprint_nts', 30)) #makes sure that if the value for how many nucleotides the ribosome cover at a time doesnt exist in the params file, then the defult will be 30 nuc.

#loading the RNA sequence files into dictionaries. 
viral_seqs = {}
for block in open(viral_fasta_path, encoding='utf-8').read().split('>')[1:]:
    head, *body = block.splitlines()
    viral_seqs[head.split('|')[0].strip()] = ''.join(body)
#calculates the total baseline dwell time for a ribosome moving through a specific RNA section
#since RNA is seq through cDNA, T nucleotides must be turnd into U.
def sequence_base_pause(virus, pos, fallback):
    """Calculate baseline ribosomal dwell time across one footprint window."""
    seq = viral_seqs.get(virus)

    if not seq:
        return fallback

    start_idx = pos - 1  # Convert 1-based biological position to 0-based Python index
    win = seq[start_idx:start_idx + FOOTPRINT_NTS].upper().replace('T', 'U')

    total = sum(
        dwell_map.get(win[i:i + 3], 0.0)
        for i in range(0, len(win) - 2, 3)
    )

    return round(total, 3) if total > 0 else fallback

#setting two constants needed for thermodynemic calculations R (gas constant) and T (the temp, human body = 37C, in kelvin)
TEMP_K = ribo_params.get("temperature_kelvin", 310.15)
R_CAL = ribo_params.get("universal_gas_constant_cal", 1.987204)

#making a dictionary for better search of the modification and and its ID
site_registry = {
    row['site_id']: {
        'subunit': row['subunit'], 'position': row['position'], 'mod_type': row['mod_type'],
        'region': row['region'], 'delta_delta_g_kcal': float(row['delta_delta_g_kcal']),
        'essentiality_alpha': float(row['essentiality_alpha']), 'data_source': row['data_source']
    }
    for _, row in df_rrna_mods.iterrows()
}

#looping through the viral data (in the future we intened to work with more viruses)
#calculating the dwell time of the ribosome on a viral "slippery" site. 
viral_registry = {}
for item in viral_metadata:
    base_pause = sequence_base_pause(
        item['virus'], 
        int(item['verified_slippery_pos_in_seq']), 
        float(item['base_pause_sec'])
    )
    viral_registry[item['virus']] = {
        'base_pause_sec': base_pause,
        'pseudoknot_unwinding_kcal': float(item['pseudoknot_unwinding_kcal']), #the energy it take to unwind the secendary structure downstream to the slippery site (adds to dwell time)
        'data_source': item['data_source']
    }

# -------------------------------------------------------------------
# 2. SITE-SPECIFIC THERAPEUTIC ACTION SPACE (EXPLICIT SITE LABELS)
# -------------------------------------------------------------------
#establishing specific modifications KO and KD to test in this pipeline (there are more than 200 modifications and looking at all possibilities of KO and KD adds too much erelevent noise to this system)
site_specific_actions = [
    {"label": "WT Control (0% Depletion)", "type": "WT", "targets": {}},
    
    # Single Site Depletions
    {"label": "28S:U1864 (PTC) 50% Depletion", "type": "Single Site Depletion", "targets": {"28S:U1864": 0.5}},
    {"label": "18S:U1248 (P-site) 50% Depletion", "type": "Single Site Depletion", "targets": {"18S:U1248": 0.5}},
    {"label": "18S:C1703 (A-site) 50% Depletion", "type": "Single Site Depletion", "targets": {"18S:C1703": 0.5}},
    {"label": "28S:U1864 (PTC) 100% Depletion", "type": "Single Site Depletion", "targets": {"28S:U1864": 1.0}},
    {"label": "18S:C1703 (A-site) 100% Depletion", "type": "Single Site Depletion", "targets": {"18S:C1703": 1.0}},
    
    # Explicit Combinatorial Site Depletions
    {"label": "Dual: U1864 + C1703 (30% Depletion)", "type": "Dual Site Depletion", "targets": {"28S:U1864": 0.3, "18S:C1703": 0.3}},
    {"label": "Dual: U1864 + C1703 (50% Depletion)", "type": "Dual Site Depletion", "targets": {"28S:U1864": 0.5, "18S:C1703": 0.5}},
    {"label": "Triple: U1864 + U1248 + C1703 (30%)", "type": "Multi Site Depletion", "targets": {"28S:U1864": 0.3, "18S:U1248": 0.3, "18S:C1703": 0.3}},
    {"label": "Pan-PTC: U1864 + U3714 + N2415 (100%)", "type": "Multi Site Depletion", "targets": {"28S:U1864": 1.0, "28S:U3714": 1.0, "28S:N2415": 1.0}}
]

# -------------------------------------------------------------------
# 3. BIOPHYSICAL KINETIC FUNCTIONS
# -------------------------------------------------------------------
#calculating the rate factor
def calculate_eyring_reduction(delta_delta_g_kcal):
    """Convert a free-energy change into a relative rate factor using the Eyring model."""
    delta_g_cal = delta_delta_g_kcal * 1000.0
    return np.exp(-delta_g_cal / (R_CAL * TEMP_K))

#simulation runs to find host safe conditions that blocks viral translation
def simulate_treatment(target_depletions, inject_noise=False):
    """Simulate host translation retention and viral frameshift-site stalling."""
    capacity_retention = 1.0 #starting point, 100% host transaltion
    #retriving factors for calculating host translation effciency
    for site_id, depletion in target_depletions.items():
        if site_id in site_registry:
            alpha = site_registry[site_id]["essentiality_alpha"]
            capacity_retention *= (1.0 - alpha * depletion)
    host_yield_pct = max(5.0, capacity_retention * 100.0) #calculating the host translation effciency

    #calculating viral inhibition
    viral_results = {}
    #Baseline Viral Features
    for virus_name, v_params in viral_registry.items():
        base_pause = v_params["base_pause_sec"]
        barrier_kcal = v_params["pseudoknot_unwinding_kcal"]
        #calculating the stall due to host "safe" modifications changes in the ribosome
        additional_stall = 0.0
        for site_id, depletion in target_depletions.items():
            if site_id in site_registry:
                ddg = site_registry[site_id]["delta_delta_g_kcal"]
                rate_factor = calculate_eyring_reduction(ddg)
                site_stall_contrib = (base_pause / rate_factor - base_pause) * (barrier_kcal / 2.0) #additional stalling due to secendary structure after the slippery site
                additional_stall += site_stall_contrib * depletion

        noise = np.random.normal(0, 0.01) if inject_noise else 0.0
        effective_pause = base_pause + additional_stall + noise #the total stalling time based on all parameteres and noise in the system
        stall_effect_pct = (additional_stall / (base_pause + additional_stall)) * 100.0 if effective_pause > base_pause else 0.0 #the total stalling that is due to the changed modifications only (no noise)
        #storing the results
        viral_results[virus_name] = {
            'slippery_pause_sec': round(max(base_pause, effective_pause), 3),
            'stall_effect_pct': round(max(0.0, stall_effect_pct), 2)
        }

    return viral_results, round(host_yield_pct, 2)

def main():
    # -------------------------------------------------------------------
    # 4. EPSILON-GREEDY CANDIDATE SEARCH & EVALUATION LOOP
    # -------------------------------------------------------------------
    #adding a fixed starting point to the random numbers
    random.seed(42)
    np.random.seed(42)

    q_table = {i: 0.0 for i in range(len(site_specific_actions))} #score board for each candidate action
    lr = 0.15 #updating weight gradually
    epsilon = 0.20 #80% of the time the next candidate will be chosen based on the highest q score, and 20% of the time it will bw random

    #running scorring simulations 800 times
    for ep in range(800):
        if random.uniform(0, 1) < epsilon: # 20% random
            action_idx = random.choice(list(q_table.keys()))
        else: # 80% highest score
            action_idx = max(q_table, key=q_table.get)
        #looking at spesific changed sites and evaluating host translation vs viral inhibition
        act = site_specific_actions[action_idx]
        viral_res, avg_host_yield = simulate_treatment(act['targets'], inject_noise=True)
        avg_viral_inhibit = np.mean([v['stall_effect_pct'] for v in viral_res.values()])
        
        host_penalty_score = (90.0 - avg_host_yield) * 20.0 if avg_host_yield < 90.0 else 0.0
        reward = avg_viral_inhibit - host_penalty_score
        
        q_table[action_idx] += lr * (reward - q_table[action_idx])

    #clean of noise evaluation for each candidate placed in a dataframe + ranking the top 5
    #this loop combains the scoring from the previuse simulation with the physical simulation metrics (exact pause durations, host yield percentages) for ranking of the candidates
    rows = []
    for idx, act in enumerate(site_specific_actions):
        viral_res, avg_host_yield = simulate_treatment(act['targets'], inject_noise=False)

        hiv_inhibit = viral_res.get('HIV-1', {}).get('stall_effect_pct', 0.0)
        sars_inhibit = viral_res.get('SARS-CoV-2', {}).get('stall_effect_pct', 0.0)
        avg_viral_inhibit = (hiv_inhibit + sars_inhibit) / 2.0 #avraging the inhibition for both viruses since they use the same mechanism
        is_safe = (avg_host_yield >= 90.0) #safe means host transaltion is kept above 90%
        #converting candidates dictionary into strings, to use for the plots
        target_str = ", ".join([f"{k} ({int(v*100)}% Depletion)" for k, v in act['targets'].items()]) if act['targets'] else "None (WT)"

        rows.append({
            'label': act['label'],
            'type': act['type'],
            'target_positions': target_str,
            'avg_host_yield_pct': avg_host_yield,
            'avg_viral_inhibit_pct': avg_viral_inhibit,
            'hiv1_inhibit_pct': hiv_inhibit,
            'sarscov2_inhibit_pct': sars_inhibit,
            'hiv1_pause_sec': viral_res.get('HIV-1', {}).get('slippery_pause_sec', 0.48),
            'sarscov2_pause_sec': viral_res.get('SARS-CoV-2', {}).get('slippery_pause_sec', 0.52),
            'q_value': round(q_table[idx], 3),
            'is_safe': is_safe
        })

    df_summary = pd.DataFrame(rows)
    top_5_safe = (
        df_summary[df_summary['is_safe']]
        .sort_values(
            by=['q_value', 'avg_host_yield_pct', 'label'],
            ascending=[False, False, True]
        )
        .head(5)
    )
    # Print Top Safe Candidates
    print("\n" + "=" * 115)
    print("        TOP HOST-SAFE SITE-SPECIFIC CANDIDATES (DATA-GROUNDED HYBRID MODEL)        ")
    print("=" * 115)
    print(top_5_safe[['label', 'target_positions', 'avg_host_yield_pct', 'avg_viral_inhibit_pct', 'q_value']].to_string(index=False))
    print("=" * 115 + "\n")

    # -------------------------------------------------------------------
    # 5. VISUALIZATION SECTION
    # -------------------------------------------------------------------
    #both plots will be presented in the same figure
    plt.figure(figsize=(18, 10))
    colors = ["#ec4899", "#3b82f6", "#8b5cf6", "#f59e0b", "#10b981"]

    # --- Panel A: Viral Inhibition vs. Host Safety Threshold ---
    ax_b = plt.subplot(1, 2, 1) #the first plot will be on the left of the grid

    # Scatter plot for all candidates
    sns.scatterplot(
        data=df_summary,
        x="avg_viral_inhibit_pct",
        y="avg_host_yield_pct",
        hue="is_safe",
        style="type",
        palette={True: "#9ca3af", False: "#ef4444"},
        s=120,
        alpha=0.5,
        ax=ax_b
    )

    # Overlay Top 5 Candidates with high-visibility markers and clean number badges
    for idx, (_, row) in enumerate(top_5_safe.iterrows()):
        x_val = row['avg_viral_inhibit_pct']
        y_val = row['avg_host_yield_pct']
        c_color = colors[idx % len(colors)]
        
        # Large colored scatter point
        ax_b.scatter(x_val, y_val, color=c_color, s=280, edgecolor="black", zorder=5)
        
        # Distinct number badge printed inside the point
        ax_b.text(x_val, y_val, f"{idx+1}", color="white", fontweight="bold", fontsize=9, ha="center", va="center", zorder=6)

    plt.axhspan(90, 102, color='#10b981', alpha=0.10)
    plt.axhline(90.0, color="#10b981", linestyle="--", alpha=0.8, linewidth=1.5)

    plt.title("Panel A: Viral Stalling Effect vs. Host Safety", fontsize=12, fontweight="bold")
    plt.xlabel("Average Viral Stalling Effect (%)", fontsize=10)
    plt.ylabel("Host Translation Yield (% of WT)", fontsize=10)

    plt.xlim([-5, 100])
    plt.ylim([70, 102])
    plt.grid(True, linestyle="--", alpha=0.5)

    # Custom Legend Handles & Labels using "Depletion"
    legend_elements = [
        Line2D([0], [0], color='#10b981', linestyle='--', lw=1.5, label='Host Safety Cutoff (≥90%)'), #candidates above that line are host "safe"
        Line2D([0], [0], marker='o', color='w', label='Safe Host (≥90%)', markerfacecolor='#9ca3af', markersize=7),
        Line2D([0], [0], marker='o', color='w', label='Toxic Host (<90%)', markerfacecolor='#ef4444', markersize=7),
        Line2D([0], [0], marker='o', color='w', label='Single Site Depletion', markerfacecolor='gray', markersize=7),
        Line2D([0], [0], marker='X', color='w', label='Dual Site Depletion', markerfacecolor='gray', markersize=7),
        Line2D([0], [0], marker='P', color='w', label='Multi Site Depletion', markerfacecolor='gray', markersize=7),
    ]

    ax_b.legend(handles=legend_elements, loc="lower left", fontsize=7.5, frameon=True, facecolor="white", framealpha=0.9)

    # Summary Table placed directly UNDER Panel A
    table_data = []
    for idx, (_, row) in enumerate(top_5_safe.iterrows()):
        table_data.append([f"#{idx+1}", row['label'], f"{row['avg_host_yield_pct']}%", f"{row['avg_viral_inhibit_pct']}%"])

    col_labels = ["Rank", "Target Action", "Host Yield", "Stalling Effect"]
    col_widths = [0.12, 0.58, 0.15, 0.15]

    tbl = plt.table(
        cellText=table_data,
        colLabels=col_labels,
        colWidths=col_widths,
        loc="bottom",
        cellLoc="center",
        bbox=[-0.05, -0.48, 1.10, 0.30]
    )
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(7.0)

    # Color the Rank column cells matching their scatter points
    for i in range(len(top_5_safe)):
        tbl[(i + 1, 0)].set_facecolor(colors[i % len(colors)])
        tbl[(i + 1, 0)].get_text().set_color("white")
        tbl[(i + 1, 0)].get_text().set_weight("bold")

    # --- Panel B: Slippery Site Stalling Duration ---
    plt.subplot(1, 2, 2) #this plot will be presented on the right of the grid

    x = np.arange(len(df_summary))
    width = 0.35

    # Fully descriptive x-axis labels corresponding to each candidate action
    short_action_labels = [
        "WT Control",
        "U1864 (50%)",
        "U1248 (50%)",
        "C1703 (50%)",
        "U1864 (100%)",
        "C1703 (100%)",
        "Dual (30%)",
        "Dual (50%)",
        "Triple (30%)",
        "Pan-PTC (100%)"
    ]

    plt.bar(x - width/2, df_summary['hiv1_pause_sec'], width, label='HIV-1 Pause Duration (s)', color='#8b5cf6')
    plt.bar(x + width/2, df_summary['sarscov2_pause_sec'], width, label='SARS-CoV-2 Pause Duration (s)', color='#f59e0b')

    plt.xticks(ticks=x, labels=short_action_labels, rotation=35, ha='right', fontsize=8)
    plt.title("Panel B: Slippery Site Stalling Duration", fontsize=12, fontweight="bold")
    plt.xlabel("Therapeutic Depletion Candidates", fontsize=10)
    plt.ylabel("Pause Duration (Seconds)", fontsize=10)
    plt.grid(True, linestyle="--", alpha=0.5)

    # Legend placed under Panel B
    plt.legend(bbox_to_anchor=(0.5, -0.22), loc="upper center", fontsize=8, frameon=True, ncol=2)

    # Adjust horizontal spacing and bottom margin
    plt.subplots_adjust(bottom=0.35, wspace=0.32)

    # Explicit definition of output path
    plot_filename = "Hybrid_MultiTarget_Annotated_Clean.png"
    plot_path = os.path.join(output_dir, plot_filename)

    # Save Plot Output
    plt.savefig(plot_path, dpi=300, bbox_inches="tight")
    plt.close()

    print(f"\n[✓] Clean plot successfully generated and saved to: {plot_path}")

if __name__ == "__main__":
    main()