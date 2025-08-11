Overview
This repository contains Stata code for the empirical analysis of the causal impact of housing prices on household fertility decisions, using microdata from the China Family Panel Studies (CFPS, 2010–2020). The scripts implement data cleaning, variable construction, descriptive statistics, baseline and instrumental-variable regressions, mechanism exploration, heterogeneity analyses, and robustness checks.

Data Sources
Household-level panel data: CFPS six waves (2010, 2012, 2014, 2016, 2018, 2020)

City-level data: Administrative statistics on water resources, population, fiscal budget, land area, schools, GDP, and credit conditions

Workflow Summary
Data Cleaning & Preparation

Standardizes household IDs and merges household, community, and city-level data

Encodes Chinese string variables (gender, health, hukou) into binary indicators

Constructs new fertility-related variables (birth indicators, age at birth, child age)

Restricts to valid birth-age range (18–40) and births from 2010 onward

Handles outliers and applies winsorization

Descriptive Statistics

Exports summary statistics for key variables to Word format

Variable Construction

Generates real income and property values (log-transformed)

Calculates community-level mean housing prices per m²

Constructs price growth measures and interaction terms for IVs

Creates city-level covariates: population, fiscal capacity, land area, water resources, schools, GDP

Econometric Analysis

Baseline Models: High-dimensional fixed effects (HDFE) regressions controlling for household, community, and city factors

Instrumental Variable Models: Instruments housing prices with fiscal capacity × water-to-land ratio

Mechanism Checks: Wealth and cost effects, homeowners vs. renters, credit constraints

Heterogeneity Analyses: By education, number of children, hukou status, child gender

Robustness Checks: Outlier removal, alternative differenced specifications

Chinese-to-English Variable Mapping
At the end of the code, a mapping dictionary links original Chinese variable names to English equivalents and definitions, ensuring transparency for replication.

Key Methods & Assumptions
Fixed Effects: Household, community, city-year FE to isolate within-unit variation

IV Relevance: City fiscal capacity and water-to-land ratio affect housing supply/cost, shifting local prices

IV Exclusion: After controlling for fixed effects and covariates, instruments affect fertility only via prices

Requirements
Stata 16+ (tested on Stata 18)

User-written packages:

winsor2 (for winsorization)

logout (for exporting tables)

reghdfe / ivreghdfe (for HDFE estimation)

Usage
Run the .do file from start to finish after placing the CFPS microdata and city-level datasets in the specified paths. Adjust use "XXXX" lines to match your local directory structure. The code will generate cleaned datasets, descriptive tables, and regression outputs.

Output
Cleaned panel datasets

Descriptive statistics (Word tables)

Regression output (HDFE and IV estimates)

Mechanism and heterogeneity results
