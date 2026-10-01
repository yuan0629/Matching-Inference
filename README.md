# Matching-Inference

## Simulation Scripts


Figure 2 (a).R examines the entrywise performance of the initial estimator obtained using Algorithm 1 and compares it with the double-enhancement algorithm of Tang et al. under one-to-one matching, across different signal magnitudes. It generates Figure 2(a).

Figure 2 (b).R conducts the same comparison of entrywise estimation errors under one-to-many matching with one-sided random arrival. It generates Figure 2(b).


Figure 2 (c).R conducts the same comparison of entrywise estimation errors under two-sided random arrival. It generates Figure 2(c).


Figure 3 (a)-(d).R evaluates the matching evaluation framework under one-to-one matching for four choices of the target matrix $Q$: a single entry, a one-to-one matching, the difference between two one-to-one matchings, and a one-to-many matching. It generates Figure 3(a)-(d).


Figure 3 (e)-(h).R evaluates the matching evaluation framework under one-to-many matching with one-sided random arrival for the same four choices of $Q$. It generates Figure 3(e)-(h).


Figure 3 (i)-(l).R evaluates the matching evaluation framework under two-sided random arrival for the same four choices of $Q$. It generates Figure 3 (i)-(l).


Appendix A.R compares the finite-sample variability and Gaussian approximation of the estimators with and without sample splitting under one-to-one matching, using $T$ in $\{1000, 2000, 5000\}$ and four choices of $Q$. It generates Table 3 and Figures 5-7 in Appendix A.


## Real data


raw data.csv was constructed using game-level umpire performance records from the [UmpScorecards data archive](https://umpscorecards.com/data/games) and game information from [the official MLB schedule](https://www.mlb.com/schedule/2025-03-27). It contains the participating teams, assigned umpires, and umpires' performance for the 2025 regular season, together with historical UmpScorecards records from 2015 to 2024.


data cleaning.R reads raw data.csv and restricts the main analysis sample to the 2025 regular season so that both sides of the matching market remain fixed throughout the study period. It retains games with an observed Accuracy Above Expected and valid team and umpire indices, verifies that each home team and each umpire appear at most once on each game day, and generates clean data.csv.



data analysis.R reads clean data.csv and applies the sample-splitting matching learning and evaluation procedure described in the paper. It generates the $p$-values for the ten dates with 15 arrived home teams shown in Table 1.
