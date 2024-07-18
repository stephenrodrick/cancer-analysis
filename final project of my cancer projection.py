import tkinter as tk
from tkinter import ttk, messagebox
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import matplotlib.pyplot as plt
from io import StringIO

# Function to analyze nucleotide sequence
def analyze_nucleotide_sequence(fasta_sequence):
    try:
        fasta_io = StringIO(fasta_sequence)
        records = list(SeqIO.parse(fasta_io, "fasta"))
        if not records:
            return "Invalid FASTA format or empty sequence."
        sequence = str(records[0].seq)
        gc_content = gc_fraction(sequence) * 100
        return f"GC Content: {gc_content:.2f}%"
    except Exception as e:
        return f"Error in analyzing nucleotide sequence: {str(e)}"

# Function to analyze protein sequence and determine cell cycle phase
def analyze_protein_sequence(protein_sequence):
    try:
        # Dummy analysis for cell cycle phase
        phases = ["G1", "S", "G2", "M"]
        import random
        cell_cycle_phase = random.choice(phases)
        phase_counts = {phase: random.randint(10, 50) for phase in phases}

        # Simulate cell cycle durations (in hours)
        phase_durations = {
            "G1": 11,
            "S": 8,
            "G2": 4,
            "M": 1
        }
        total_duration = sum(phase_durations.values())
        return cell_cycle_phase, phase_counts, phase_durations, total_duration
    except Exception as e:
        return f"Error in analyzing protein sequence: {str(e)}"

# Function to display pie chart of cell cycle phases
def display_pie_chart(phase_counts, phase_durations, total_duration):
    try:
        phases = list(phase_counts.keys())
        counts = list(phase_counts.values())
        durations = [phase_durations[phase] for phase in phases]
        
        fig, ax = plt.subplots(1, 2, figsize=(12, 6))
        
        # Pie chart for phase distribution
        ax[0].pie(counts, labels=phases, autopct='%1.1f%%')
        ax[0].set_title('Cell Cycle Phases Distribution')
        
        # Bar chart for phase durations
        ax[1].bar(phases, durations)
        ax[1].set_xlabel('Cell Cycle Phase')
        ax[1].set_ylabel('Duration (hours)')
        ax[1].set_title(f'Cell Cycle Phase Durations (Total: {total_duration} hours)')
        
        plt.show()
    except Exception as e:
        messagebox.showerror("Error", f"Error in displaying charts: {str(e)}")

# Function to handle Analyze button click
def analyze_sequences():
    nucleotide_sequence = nucleotide_text.get("1.0", tk.END).strip()
    protein_sequence = protein_text.get("1.0", tk.END).strip()

    if nucleotide_sequence:
        nucleotide_result = analyze_nucleotide_sequence(nucleotide_sequence)
        nucleotide_result_label.config(text=nucleotide_result)
    else:
        messagebox.showwarning("Input Error", "Please enter a nucleotide sequence.")

    if protein_sequence:
        result = analyze_protein_sequence(protein_sequence)
        if isinstance(result, str):
            protein_result_label.config(text=result)
        else:
            protein_result, phase_counts, phase_durations, total_duration = result
            protein_result_label.config(text=f"Cell Cycle Phase: {protein_result}")
            display_pie_chart(phase_counts, phase_durations, total_duration)
    else:
        messagebox.showwarning("Input Error", "Please enter a protein sequence.")

# Setting up the GUI
root = tk.Tk()
root.title("Cancer Cell Analysis")

# Frame for nucleotide sequence input
nucleotide_frame = ttk.LabelFrame(root, text="Nucleotide FASTA Sequence Input")
nucleotide_frame.grid(row=0, column=0, padx=10, pady=10, sticky="ew")

nucleotide_text = tk.Text(nucleotide_frame, height=10, width=50)
nucleotide_text.grid(row=0, column=0, padx=5, pady=5)

# Frame for protein sequence input
protein_frame = ttk.LabelFrame(root, text="Protein Sequence Input")
protein_frame.grid(row=1, column=0, padx=10, pady=10, sticky="ew")

protein_text = tk.Text(protein_frame, height=5, width=50)
protein_text.grid(row=0, column=0, padx=5, pady=5)

# Analyze button
analyze_button = ttk.Button(root, text="Analyze", command=analyze_sequences)
analyze_button.grid(row=2, column=0, padx=10, pady=10)

# Results display
results_frame = ttk.LabelFrame(root, text="Results")
results_frame.grid(row=3, column=0, padx=10, pady=10, sticky="ew")

nucleotide_result_label = ttk.Label(results_frame, text="Nucleotide Sequence Analysis Result")
nucleotide_result_label.grid(row=0, column=0, padx=5, pady=5)

protein_result_label = ttk.Label(results_frame, text="Protein Sequence Analysis Result")
protein_result_label.grid(row=1, column=0, padx=5, pady=5)

# Start the GUI event loop
root.mainloop()
