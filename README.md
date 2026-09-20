# Athlete-Optimized VO2 Prediction from Minimal Physiological Signals

This project develops and evaluates models for predicting oxygen consumption and VO2max in athletes using minimal physiological and exercise-related signals.

The main objective is to estimate the time-varying relative oxygen consumption of an athlete from heart rate, running speed, and individual characteristics. The models are validated against laboratory measurements obtained with a COSMED K5 metabolic analyzer.

## Project Overview

The project is based on an incremental treadmill exercise protocol including:

- A warm-up phase.
- Progressive speed stages.
- A maximal or exhaustion phase.
- A post-exercise recovery phase.

The analysis includes several complementary approaches:

- Prediction of the VO2 trajectory over time.
- Mixed-effects modeling for repeated measurements from multiple athletes.
- Leave-one-athlete-out validation for unseen-athlete prediction.
- VO2max estimation from energy expenditure during the final high-intensity phase.
- VO2max estimation from heart-rate recovery.
- Analysis of heart-rate behavior and cardiovascular economy.
- Comparison with standard MET and ACSM-based equations.

## Main Modeling Approach

The primary model predicts relative oxygen consumption using:

- Relative running speed.
- Relative heart rate.
- Age.
- Body weight.
- Height.
- Athlete-specific variability through mixed-effects modeling.

The fixed-effects component can be applied to a new athlete when individual random effects are not yet available. Because heart rate and relative speed vary over time, the model can estimate a complete VO2 trajectory rather than only a single value.

## Validation

The model is evaluated using leave-one-athlete-out cross-validation. In this procedure, one athlete is excluded from model fitting and then used as an unseen test subject.

This validation strategy is designed to assess how well the model generalizes to athletes who were not included during training.

## Requirements

The analysis requires:

- MATLAB.
- Statistics and Machine Learning Toolbox.
- Input exercise-test files in `.xlsx`, `.xls`, or `.csv` format.
- Laboratory VO2 measurements and heart-rate data organized according to the project data format.

## Usage

Open MATLAB and run:

```matlab
vo2_mixed_effects_analysis
```

The script opens a file-selection window. Select one or more exercise-test data files and follow the instructions shown in the MATLAB Command Window.

The script processes the selected tests, extracts the exercise stages, calculates physiological features, fits the models, performs validation, generates diagnostic figures, and saves the analysis results as:

```text
vo2_analysis_results.xlsx
```

## Data and Privacy

Raw athlete data, personal information, laboratory files, and generated Excel results are not included in this repository. These files may contain sensitive physiological and anthropometric information and should be stored locally or in an appropriately protected location.

## Scope

This repository contains research code for an athlete-focused VO2 and VO2max prediction workflow. The models are intended for research and experimental analysis and should not be interpreted as clinical diagnostic tools.
