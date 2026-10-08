# rRNA Modification and Viral Translation Model

### Biophysical Simulation and Epsilon-Greedy Screening of rRNA Modification Depletions

## Project Overview

This project started with a simple goal — to predict the outcome of a biological experiment before performing it in the lab.

I wanted to investigate how changes in ribosomal RNA (rRNA) modifications might influence viral protein translation. However, I couldn't find an existing computational model that directly incorporated the specific rRNA modifications I wanted to investigate.

This led me to develop a Python-based biophysical simulation framework exploring the effects of rRNA modification depletion on ribosomal pausing at programmed frameshifting sites in HIV-1 and SARS-CoV-2.

The model also considers how these changes might affect normal host translation.

## How the Model Works

The framework combines several computational approaches.

- **Biophysical simulation** using thermodynamic assumptions and ribosomal dwell-time inputs to estimate changes in ribosomal pausing.
- **rRNA modification screening** evaluating selected single-site and combined modification-depletion scenarios.
- **Epsilon-greedy optimization** ranking candidate scenarios according to their predicted viral stalling effects and modeled host translation retention.
- **Data visualization** comparing predicted ribosomal behavior across different conditions.

The model evaluates ten predefined scenarios, including an unmodified wild-type control.

## Input Data and Limitations

The current model uses a combination of biological reference information, physical constants, estimated parameters, and simulated values.

The input files include

- `riboseq_dwell_times.csv` — Codon-specific dwell-time inputs.
- `ribosome_params.json` — Ribosomal and thermodynamic parameters.
- `rrna_modifications_curated.csv` — Selected rRNA modification sites and modeled energetic and functional parameters.
- `viral_frameshift_sites.json` — Viral frameshifting sites and associated model parameters.
- `viral_sequences.fasta` — Viral nucleotide sequences.

**Important** — This is an exploratory computational model. Several input parameters are assumed or estimated rather than experimentally measured, and their source annotations require verification.

The 90% host translation retention threshold is a modeling choice, not an experimentally established safety threshold.

The current predictions have not been experimentally validated and should not be interpreted as demonstrated viral inhibition or measured changes in frameshifting efficiency.

## Installation and Usage

The project was developed using Python 3.12.

Install the required libraries.

```bash
pip install numpy pandas matplotlib seaborn
```

Keep `main.py` and all five input files in the same directory.

Run the program.

```bash
python main.py
```

The program prints a table of the five highest-ranking candidates that retain at least 90% of modeled host translation capacity.

It also generates a two-panel figure comparing

- **Panel A** — Average modeled viral stalling effect versus host translation retention.
- **Panel B** — Predicted ribosomal pause duration at HIV-1 and SARS-CoV-2 frameshifting sites.

The figure is saved to `output/Hybrid_MultiTarget_Annotated_Clean.png`.

## Reproducibility

Python and NumPy random seeds are fixed to 42 to support reproducibility of the candidate-ranking procedure.

## Future Development

The next stage of this project will focus on incorporating real experimental measurements to replace or calibrate estimated parameters and evaluate the model's predictive accuracy.

I also plan to explore more advanced machine learning approaches to better characterize the relationships between rRNA modifications, ribosomal dynamics, and viral translation.

The long-term goal is to develop a computational framework that can generate testable biological predictions and help guide future laboratory experiments.

## Project Status

**Under active development — experimental validation planned.**
