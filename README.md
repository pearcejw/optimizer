# Portfolio Optimizer

A strategic asset allocation tool built in R Shiny. The app combines mean-variance optimization, Black-Litterman, and Monte Carlo simulation into a single interactive dashboard for building and analysing multi-asset portfolios.

------------------------------------------------------------------------

## Requirements

``` r
install.packages(c(
  "shiny", "shinydashboard", "DT", "quadprog",
  "plotly", "dplyr", "tidyr", "shinyjs", "MASS", "writexl"
))
```

R 4.1 or later is recommended.

------------------------------------------------------------------------

## Running the App

``` r
shiny::runApp("app.R")
```

Or open `app.R` in RStudio and click **Run App**.

------------------------------------------------------------------------

## App Structure

The app is organised into nine sidebar tabs, intended to be worked through roughly in order.

### 1. Assets & Assumptions

The starting point. Define the universe of assets and their capital market assumptions.

-   **Risk-Free Rate** — sets the rate used for Sharpe ratio calculations throughout the app
-   **Asset table** — edit expected returns (% p.a.), volatility (% p.a.), and asset class directly in the table; also includes a **Current (%)** column for entering your existing portfolio weights
-   **Add / Remove assets** — extend or trim the universe; the correlation matrix adjusts automatically
-   **Paste from spreadsheet** — paste directly from Excel or Google Sheets; columns are auto-detected by name (Asset, Return, Vol, Class, Current); tab- and comma-separated formats both work

### 2. Correlations

-   Edit the pairwise correlation matrix directly in the table (values −100 to +100)
-   The matrix is kept symmetric automatically; the diagonal is locked at 100
-   **Paste from spreadsheet** — paste a labelled square matrix; rows/columns are reordered to match the asset list if names are recognisable

### 3. Constraints

-   **Global weight bounds** — set minimum and maximum weight per asset globally; negative minimums allow short positions
-   **Per-asset overrides** — fine-tune bounds for individual assets; blank fields inherit the global setting
-   **Asset class constraints** — cap or floor the total allocation to an entire class (e.g. all Equity between 30–60%)
-   **Custom group constraints** — combine any subset of assets into a group with its own allocation bounds
-   **Solver settings** — control the number of points computed on the efficient frontier

### 4. Efficient Frontier

Click **Run Optimizer** to compute the frontier. The chart plots the full frontier, individual assets, the Global Minimum Variance (GMV) portfolio, and the Maximum Sharpe Ratio (MSR) portfolio. If a Current Portfolio has been entered it appears as an orange marker, allowing direct visual comparison against the frontier.

Summary statistics (return, volatility, Sharpe) are shown for the GMV, MSR, and any portfolio selected in the Analysis tab.

### 5. Portfolio Analysis

Select an objective function and click **Analyse** to evaluate a specific portfolio. The objective function dropdown includes:

| Objective | Description |
|----|----|
| Minimum Variance | Lowest possible volatility |
| Maximum Sharpe Ratio | Highest risk-adjusted return |
| Maximum Diversification | Maximizes the diversification ratio |
| Risk Parity | Equal risk contribution from each asset |
| Equal Weight | Naïve 1/N allocation |
| Min CVaR (5%) | Minimizes expected shortfall at the 5% level |
| Maximum Return | Highest achievable return given constraints |
| Maximum Utility | Mean-variance utility with configurable risk aversion λ |
| Min Tracking Error | Minimizes tracking error against a user-defined benchmark |
| Target Return | Min-variance portfolio for a specified return target |
| Target Volatility | Portfolio closest to a specified volatility target |
| Manual Weights | Enter weights directly; normalized to sum to 100% automatically |
| Current Portfolio | Analyses the portfolio entered in the Assets tab |

Results include value boxes for return, volatility, Sharpe, and CVaR; allocation and risk contribution bar charts; and a detailed weights table.

### 6. Portfolio Comparison

Compare up to three portfolios side by side. Each slot (A, B, C) has its own objective dropdown drawing from the same list as the Analysis tab plus a "Current Analysis" option that mirrors whatever was last run on the Analysis tab. Click **Compare** to generate:

-   **Stats cards** — return, volatility, Sharpe, and CVaR for each portfolio, colour-coded orange / blue / green
-   **Grouped bar chart** — weights across all assets for all selected portfolios
-   **Side-by-side table** — weight and risk contribution columns for each portfolio

### 7. Black-Litterman

Blends CAPM-implied equilibrium returns with subjective investor views using Bayesian updating. The posterior return vector replaces raw expected returns in the optimizer.

-   **Market-cap weights** — sets the prior equilibrium
-   **Risk aversion (δ)** and **prior confidence (τ)** — control the strength of the equilibrium prior
-   **Views** — each view specifies a P-vector (long/short asset weights), an expected return Q, and a confidence level (0 = ignore, 1 = certain)

Click **Compute BL Posterior** to see the equilibrium vs. posterior return comparison and the BL Max-Sharpe portfolio allocation.

### 8. Monte Carlo

Simulates portfolio wealth paths using a multivariate normal distribution over a multi-year horizon.

**Settings:** portfolio to simulate, horizon, number of simulations, initial value, annual contribution, and CVaR confidence level.

**Outputs:**

-   **Wealth path chart** — fan of individual paths with 5th/25th/50th/75th/95th percentile bands
-   **Value boxes** — median terminal value, 10th and 90th percentile terminal values, expected loss at the CVaR level
-   **Terminal Wealth Distribution** — histogram of final portfolio values
-   **Annual Return Boxplots** — distribution of single-year returns by year
-   **Annualised Return Distribution (CAGR)** — kernel density of full-horizon compound annualised returns across all paths, with shaded IQR and tail regions and labelled percentile markers
-   **Tail Risk Summary** — table showing best/worst annualised return and worst-case maximum drawdown at the 90%, 95%, and 99% confidence levels

### 9. Comparison & Export

Generates a full comparison table across all systematic portfolios (Min Variance, Max Sharpe, Max Diversification, Risk Parity, Equal Weight, Max Utility, Min CVaR, Max Return) plus the Current Portfolio, Black-Litterman, and Current Analysis portfolios if available.

Statistics shown: expected return, volatility, Sharpe ratio, CVaR (5%), and diversification ratio.

A weights heatmap visualises allocation differences across all portfolios.

**Export options:**

| Button | Contents |
|----|----|
| Export CSV | Weights + statistics for all portfolios |
| Export Excel | Multi-sheet workbook: Weights, Statistics, Assets, Correlations, Black-Litterman (if run) |
| Export BL CSV | Black-Litterman equilibrium vs. posterior returns and weights |
| Export MC Paths | Full matrix of simulated wealth paths (one row per simulation) |

------------------------------------------------------------------------

## Current Portfolio Feature

The **Current (%)** column in the Assets table lets you record your existing portfolio allocation. Weights do not need to sum to 100 — the app normalises them automatically. Dollar values or relative proportions work equally well.

Once entered, the current portfolio:

-   Appears as an **orange square marker** on the Efficient Frontier chart
-   Is available as an objective in the **Portfolio Analysis** and **Portfolio Comparison** tabs
-   Pre-populates the **Min Tracking Error** benchmark inputs
-   Is included as a row in the **Comparison & Export** table
-   Can be selected for **Monte Carlo** simulation

------------------------------------------------------------------------

## Notes on the Optimizer

-   All optimization uses quadratic programming via the `quadprog` package
-   The covariance matrix is built from volatilities and the correlation matrix; no historical return series are required
-   Constraints are passed to every optimization call — the frontier, key portfolios, and all analysis results all respect the same bounds
-   Risk Parity uses an iterative equal-risk-contribution algorithm; it does not respect weight constraints by design
-   Min CVaR uses simulation (5,000 draws) over the efficient frontier to find the minimum-CVaR portfolio; it is the slowest objective to compute

------------------------------------------------------------------------

## File Structure

```         
optimizer/
├── app.R          # Single-file Shiny application (UI + server)
├── README.md      # This file
└── optimizer.Rproj
```
