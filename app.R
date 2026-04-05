# ============================================================
# Strategic Asset Allocation Optimizer — Extended Edition
# ============================================================
# Required packages:
#   shiny, shinydashboard, DT, quadprog, plotly,
#   dplyr, tidyr, shinyjs, MASS, writexl
#
# install.packages(c("shiny","shinydashboard","DT","quadprog",
#   "plotly","dplyr","tidyr","shinyjs","MASS","writexl"))
# ============================================================

# Load Required Libraries -------------------------------------------------

library(shiny)
library(shinydashboard)
library(DT)
library(quadprog)
library(plotly)
library(dplyr)
library(tidyr)
library(shinyjs)
library(MASS)
library(writexl)

# Define Helper Functions -------------------------------------------------

## If Else Function --------------------------------------------------------

`%||%` <- function(a, b) if (!is.null(a)) a else b

## Build Covariance Matrix -------------------------------------------------

build_cov_matrix <- function(vols, cor_mat) {
  D <- diag(vols / 100)
  C <- cor_mat / 100; diag(C) <- 1
  D %*% C %*% D
}

## QP Solver ---------------------------------------------------------------
# Min-variance with optional return target

solve_mv_portfolio <- function(mu, Sigma, target_return = NULL,
                               w_min_vec = NULL, w_max_vec = NULL,
                               group_constraints = NULL) {
  n <- length(mu)
  if (is.null(w_min_vec)) w_min_vec <- rep(0, n)
  if (is.null(w_max_vec)) w_max_vec <- rep(1, n)
  Dmat <- 2 * Sigma + diag(1e-8, n)
  dvec <- rep(0, n)
  Aeq  <- matrix(1, 1, n); beq <- 1
  if (!is.null(target_return)) { Aeq <- rbind(Aeq, mu); beq <- c(beq, target_return) }
  Amat <- t(rbind(Aeq, diag(n), -diag(n)))
  bvec <- c(beq, w_min_vec, -w_max_vec)
  meq  <- nrow(Aeq)
  if (!is.null(group_constraints)) {
    for (gc in group_constraints) {
      if (is.null(gc) || length(gc$indices) == 0) next
      rm <- rep(0, n); rm[gc$indices] <-  1
      rx <- rep(0, n); rx[gc$indices] <- -1
      Amat <- cbind(Amat, rm, rx)
      bvec <- c(bvec, gc$min/100, -gc$max/100)
    }
  }
  tryCatch({
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
    w   <- sol$solution; w[abs(w) < 1e-7] <- 0; w / sum(w)
  }, error = function(e) NULL)
}

## Efficient Frontier ------------------------------------------------------


generate_frontier <- function(mu, Sigma, n_points = 60,
                              w_min_vec = NULL, w_max_vec = NULL,
                              group_constraints = NULL) {
  n <- length(mu)
  if (is.null(w_min_vec)) w_min_vec <- rep(0, n)
  if (is.null(w_max_vec)) w_max_vec <- rep(1, n)
  w_gmv <- solve_mv_portfolio(mu, Sigma, w_min_vec = w_min_vec,
                              w_max_vec = w_max_vec, group_constraints = group_constraints)
  if (is.null(w_gmv)) return(NULL)
  mu_min <- sum(w_gmv * mu)
  mu_max <- min(max(mu), sum(w_max_vec * sort(mu, decreasing = TRUE)))
  if (mu_max <= mu_min) mu_max <- mu_min + 0.02
  results <- lapply(seq(mu_min, mu_max, length.out = n_points), function(r) {
    w <- solve_mv_portfolio(mu, Sigma, target_return = r,
                            w_min_vec = w_min_vec, w_max_vec = w_max_vec,
                            group_constraints = group_constraints)
    if (is.null(w)) return(NULL)
    list(Return = sum(w * mu),
         Vol    = sqrt(as.numeric(t(w) %*% Sigma %*% w)) * 100,
         Weights = w)
  })
  results[!sapply(results, is.null)]
}

## Objective Functions -----------------------------------------------------

max_sharpe_portfolio <- function(frontier, rf = 0) {
  sh <- sapply(frontier, function(p) (p$Return - rf) / (p$Vol / 100))
  frontier[[which.max(sh)]]
}

max_div_portfolio <- function(Sigma, vols, w_min_vec = NULL, w_max_vec = NULL,
                              group_constraints = NULL) {
  solve_mv_portfolio(vols / 100, Sigma, w_min_vec = w_min_vec,
                     w_max_vec = w_max_vec, group_constraints = group_constraints)
}

risk_parity_portfolio <- function(Sigma, max_iter = 2000, tol = 1e-10) {
  n <- nrow(Sigma); w <- rep(1/n, n)
  for (i in seq_len(max_iter)) {
    sp <- sqrt(as.numeric(t(w) %*% Sigma %*% w)) # Calculate portfolio vol
    RC <- w * (as.numeric(Sigma %*% w) / sp) # Calculate contribution to vol
    # wn <- w / (n * RC / sp); wn <- wn / sum(wn) # Original risk-parity formula
    wn <- w / RC; wn <- wn / sum(wn) # Update to fix risk-parity formula
    if (max(abs(wn - w)) < tol) { w <- wn; break }
    w <- wn
  }
  w
}

equal_weight_portfolio <- function(n) rep(1/n, n)

min_cvar_portfolio <- function(mu, Sigma, alpha = 0.05, n_sim = 5000,
                               w_min_vec = NULL, w_max_vec = NULL,
                               group_constraints = NULL) {
  set.seed(42)
  sc   <- MASS::mvrnorm(n_sim, mu = mu, Sigma = Sigma)
  front <- generate_frontier(mu, Sigma, n_points = 80,
                             w_min_vec = w_min_vec, w_max_vec = w_max_vec,
                             group_constraints = group_constraints)
  if (is.null(front)) return(NULL)
  cv <- sapply(front, function(p) {
    pr <- as.numeric(sc %*% p$Weights)
    -mean(pr[pr <= quantile(pr, alpha)])
  })
  front[[which.min(cv)]]$Weights
}

max_return_portfolio <- function(mu, w_min_vec = NULL, w_max_vec = NULL) {
  n <- length(mu)
  if (is.null(w_min_vec)) w_min_vec <- rep(0, n)
  if (is.null(w_max_vec)) w_max_vec <- rep(1, n)
  w <- w_min_vec; rem <- 1 - sum(w)
  for (i in order(mu, decreasing = TRUE)) {
    add <- min(rem, w_max_vec[i] - w[i]); w[i] <- w[i] + add; rem <- rem - add
    if (rem < 1e-9) break
  }
  w / sum(w)
}

max_utility_portfolio <- function(mu, Sigma, lambda = 3,
                                  w_min_vec = NULL, w_max_vec = NULL,
                                  group_constraints = NULL) {
  n <- length(mu)
  if (is.null(w_min_vec)) w_min_vec <- rep(0, n)
  if (is.null(w_max_vec)) w_max_vec <- rep(1, n)
  Dmat <- lambda * Sigma + diag(1e-8, n); dvec <- mu
  Aeq <- matrix(1, 1, n); beq <- 1
  Amat <- t(rbind(Aeq, diag(n), -diag(n)))
  bvec <- c(beq, w_min_vec, -w_max_vec); meq <- 1
  if (!is.null(group_constraints)) {
    for (gc in group_constraints) {
      if (is.null(gc) || length(gc$indices) == 0) next
      rm <- rep(0, n); rm[gc$indices] <-  1
      rx <- rep(0, n); rx[gc$indices] <- -1
      Amat <- cbind(Amat, rm, rx); bvec <- c(bvec, gc$min/100, -gc$max/100)
    }
  }
  tryCatch({
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
    w   <- sol$solution; w[abs(w) < 1e-7] <- 0; w / sum(w)
  }, error = function(e) NULL)
}

min_te_portfolio <- function(mu, Sigma, w_bench, w_min_vec = NULL,
                             w_max_vec = NULL, group_constraints = NULL) {
  n <- length(mu)
  if (is.null(w_min_vec)) w_min_vec <- rep(0, n)
  if (is.null(w_max_vec)) w_max_vec <- rep(1, n)
  Dmat <- 2 * Sigma + diag(1e-8, n)
  dvec <- 2 * as.numeric(Sigma %*% w_bench)
  Aeq  <- matrix(1, 1, n); beq <- 1
  Amat <- t(rbind(Aeq, diag(n), -diag(n)))
  bvec <- c(beq, w_min_vec, -w_max_vec); meq <- 1
  if (!is.null(group_constraints)) {
    for (gc in group_constraints) {
      if (is.null(gc) || length(gc$indices) == 0) next
      rm <- rep(0, n); rm[gc$indices] <-  1
      rx <- rep(0, n); rx[gc$indices] <- -1
      Amat <- cbind(Amat, rm, rx); bvec <- c(bvec, gc$min/100, -gc$max/100)
    }
  }
  tryCatch({
    sol <- solve.QP(Dmat, dvec, Amat, bvec, meq = meq)
    w   <- sol$solution; w[abs(w) < 1e-7] <- 0; w / sum(w)
  }, error = function(e) NULL)
}

## ── Black-Litterman ───────────────────────────────────────────
bl_implied_returns <- function(Sigma, w_mkt, delta = 2.5, rf = 0)
  as.numeric(delta * Sigma %*% w_mkt) + rf

bl_posterior <- function(mu_eq, Sigma, P, Q, omega = NULL, tau = 0.05) {
  if (!is.matrix(P)) P <- matrix(P, nrow = 1)
  if (is.null(omega)) omega <- tau * P %*% Sigma %*% t(P)
  tS <- tau * Sigma
  M1 <- solve(tS) + t(P) %*% solve(omega) %*% P
  M2 <- solve(tS) %*% mu_eq + t(P) %*% solve(omega) %*% Q
  list(mu = as.numeric(solve(M1, M2)), Sigma = Sigma + solve(M1))
}

## Portfolio Statistics ----------------------------------------------------

port_stats <- function(w, mu, Sigma, rf = 0) {
  ret <- sum(w * mu)
  vol <- sqrt(as.numeric(t(w) %*% Sigma %*% w)) * 100
  list(Return = ret, Vol = vol,
       Sharpe = if (vol > 0) (ret - rf) / (vol / 100) else NA_real_)
}

risk_contributions <- function(w, Sigma) {
  sp <- sqrt(as.numeric(t(w) %*% Sigma %*% w))
  RC <- w * (as.numeric(Sigma %*% w) / sp)
  # list(RC_pct = RC / sum(abs(RC)) * 100, vol = sp * 100)
  list(RC_pct = RC / sum(RC) * 100, vol = sp * 100) # Recommended change
}

## Monte Carlo -------------------------------------------------------------

run_mc_sim <- function(w, mu, Sigma, horizon, n_sim, initial, contrib = 0) {
  set.seed(123)
  raw <- MASS::mvrnorm(n_sim * horizon, mu = mu, Sigma = Sigma)
  pr  <- matrix(as.numeric(raw %*% w), nrow = n_sim, ncol = horizon)
  wl  <- matrix(0, nrow = n_sim, ncol = horizon + 1)
  wl[, 1] <- initial
  for (yr in seq_len(horizon)) wl[, yr+1] <- wl[, yr] * (1 + pr[, yr]) + contrib
  list(wealth = wl, port_ret = pr)
}

fmt_dollar <- function(x) paste0("$", format(round(x), big.mark = ",", scientific = FALSE))

# Default Portfolios ------------------------------------------------------

default_assets <- data.frame(
  Asset       = c("US Equity","Intl Equity","EM Equity",
                  "US Bonds","Intl Bonds","Real Assets","Cash"),
  Return_pct  = c(7.5, 6.5, 8.0, 3.0, 2.5, 5.0, 1.5),
  Vol_pct     = c(16.0, 18.0, 22.0, 5.0, 7.0, 14.0, 0.5),
  Class       = c("Equity","Equity","Equity",
                  "Fixed Income","Fixed Income","Alternatives","Cash"),
  Current_pct = c(30, 15, 5, 20, 15, 14, 1),
  stringsAsFactors = FALSE
)
default_cor <- matrix(c(
  100, 75, 65,-10,  5, 40,  0,
  75,100, 70,-15,  0, 35,  0,
  65, 70,100,-20, -5, 30,  0,
  -10,-15,-20,100, 60, 10, 10,
  5,  0, -5, 60,100, 10,  5,
  40, 35, 30, 10, 10,100,  0,
  0,  0,  0, 10,  5,  0,100
), nrow = 7, byrow = TRUE)

PAL <- c("#c9a84c","#4a90d9","#2ecc71","#e74c3c","#9b59b6",
         "#1abc9c","#e67e22","#95a5a6","#f1c40f","#3498db",
         "#e91e8c","#00bcd4","#8bc34a","#ff5722","#607d8b")

# CSS ---------------------------------------------------------------------

# Palette aligned with Bootswatch Cosmo (dark-mode adaptation):
#   --bg:        #0f1117   page background
#   --surface:   #1a1f2e   card / sidebar surfaces  (Cosmo navy darkened)
#   --surface2:  #222840   elevated surface / header rows
#   --border:    #2d3555   border / divider          (Cosmo border tone)
#   --primary:   #2780e3   Cosmo $primary
#   --accent:    #3498db   secondary interactive
#   --text:      #d0d6e8   body text                 (Cosmo light text)
#   --muted:     #7a87a8   muted / label text
#   --mono:      'IBM Plex Mono', 'Courier New', monospace
#   --sans:      'Source Sans 3', 'Segoe UI', system-ui, sans-serif
#                          (Cosmo uses Source Sans; IBM Plex kept as mono)
APP_CSS <- "
@import url('https://fonts.googleapis.com/css2?family=Source+Sans+3:wght@300;400;600&family=IBM+Plex+Mono:wght@400;500&display=swap');

/* ── Custom properties (single source of truth) ─────────── */
:root {
  --bg:       #0f1117;
  --surface:  #1a1f2e;
  --surface2: #222840;
  --border:   #2d3555;
  --primary:  #2780e3;
  --accent:   #3d9be9;
  --text:     #d0d6e8;
  --muted:    #7a87a8;
  --success:  #3fb618;
  --danger:   #ff0039;
  --warning:  #c9a84c;
  --sans:     'Source Sans 3', 'Segoe UI', system-ui, sans-serif;
  --mono:     'IBM Plex Mono', 'Courier New', monospace;
  --radius:   4px;
}

/* ── Base ───────────────────────────────────────────────── */
body, .content-wrapper, .main-footer {
  background: var(--bg) !important;
  font-family: var(--sans);
  color: var(--text);
}

/* ── AdminLTE chrome ────────────────────────────────────── */
.skin-blue .main-header .logo,
.skin-blue .main-header .navbar     { background: var(--surface) !important; border-bottom: 1px solid var(--border); }
.skin-blue .main-sidebar            { background: var(--surface) !important; border-right:  1px solid var(--border); }
.skin-blue .sidebar-menu > li > a   { color: var(--muted) !important; font-size: 13px; }
.skin-blue .sidebar-menu > li.active > a,
.skin-blue .sidebar-menu > li > a:hover {
  color: var(--text) !important;
  background: var(--surface2) !important;
  border-left: 3px solid var(--primary);
}

/* ── Boxes ──────────────────────────────────────────────── */
.box        { background: var(--surface) !important; border: 1px solid var(--border) !important; border-radius: var(--radius) !important; box-shadow: none !important; }
.box-header { background: var(--surface2) !important; border-bottom: 1px solid var(--border) !important; border-radius: var(--radius) var(--radius) 0 0 !important; padding: 10px 15px !important; }
.box-title  { color: var(--text) !important; font-family: var(--sans) !important; font-size: 14px !important; font-weight: 600; letter-spacing: .2px; }

/* ── Forms ──────────────────────────────────────────────── */
label                       { color: var(--muted) !important; font-size: 12px; }
.form-control               { background: var(--bg) !important; border: 1px solid var(--border) !important; color: var(--text) !important; border-radius: var(--radius) !important; font-size: 13px; }
.form-control:focus         { border-color: var(--primary) !important; box-shadow: 0 0 0 2px rgba(39,128,227,.2) !important; }
select.form-control option  { background: var(--surface); color: var(--text); }

/* ── Buttons ────────────────────────────────────────────── */
.btn-primary        { background: var(--primary) !important; border-color: var(--primary) !important; color: #fff !important; font-size: 13px; font-weight: 600; border-radius: var(--radius); }
.btn-primary:hover  { background: var(--accent) !important; border-color: var(--accent) !important; }
.btn-default        { background: var(--surface2) !important; border: 1px solid var(--border) !important; color: var(--muted) !important; font-size: 12px; border-radius: var(--radius); }
.btn-default:hover  { border-color: var(--primary) !important; color: var(--text) !important; }
.btn-danger         { background: var(--surface2) !important; border: 1px solid var(--danger) !important; color: var(--danger) !important; font-size: 12px; border-radius: var(--radius); }

/* ── DataTables ─────────────────────────────────────────── */
.dataTables_wrapper                               { color: var(--muted); font-family: var(--mono); font-size: 12px; }
.dataTables_wrapper .dataTables_length select,
.dataTables_wrapper .dataTables_filter input      { background: var(--bg) !important; color: var(--text) !important; border: 1px solid var(--border) !important; }
table.dataTable                                   { background: var(--surface) !important; color: var(--text) !important; border: none !important; }
table.dataTable thead th                          { background: var(--surface2) !important; color: var(--text) !important; border-bottom: 1px solid var(--border) !important; font-family: var(--sans); font-size: 11px; text-transform: uppercase; letter-spacing: .7px; font-weight: 600; }
table.dataTable tbody tr                          { border-bottom: 1px solid var(--surface2) !important; }
table.dataTable tbody tr:hover                    { background: var(--surface2) !important; }
table.dataTable input                             { background: var(--bg); border: 1px solid var(--border); color: var(--text); font-family: var(--mono); font-size: 12px; border-radius: var(--radius); padding: 2px 6px; width: 100%; box-sizing: border-box; }
table.dataTable input:focus                       { border-color: var(--primary); outline: none; }

/* ── Tabs ───────────────────────────────────────────────── */
.nav-tabs                  { border-bottom: 1px solid var(--border) !important; }
.nav-tabs > li > a         { color: var(--muted) !important; font-size: 13px; background: transparent !important; border: none !important; border-bottom: 2px solid transparent !important; margin-bottom: -1px; padding: 8px 15px; }
.nav-tabs > li.active > a  { color: var(--text) !important; border-bottom: 2px solid var(--primary) !important; background: transparent !important; }
.nav-tabs > li > a:hover   { color: var(--text) !important; }

/* ── Value boxes ────────────────────────────────────────── */
.value-box                      { border-radius: var(--radius) !important; border: 1px solid var(--border) !important; }
.value-box .value-box-value     { font-family: var(--mono) !important; font-size: 22px !important; color: var(--text) !important; }
.value-box .value-box-description { font-size: 11px !important; color: var(--muted) !important; text-transform: uppercase; letter-spacing: .4px; }

/* ── Sliders (ionRangeSlider) ───────────────────────────── */
.irs-bar, .irs-bar-edge  { background: var(--primary) !important; border: none !important; }
.irs-single              { background: var(--primary) !important; color: #fff !important; font-family: var(--mono) !important; font-size: 11px; }
.irs-min, .irs-max       { color: var(--muted) !important; font-size: 10px; }

/* ── Wells ──────────────────────────────────────────────── */
.well     { background: var(--surface2) !important; border: 1px solid var(--border) !important; border-radius: var(--radius) !important; padding: 12px !important; margin-bottom: 10px !important; }
.well-sm  { padding: 8px !important; }

/* ── Status badges ──────────────────────────────────────── */
.status-badge  { display: inline-block; padding: 2px 8px; border-radius: var(--radius); font-family: var(--mono); font-size: 11px; }
.status-ok     { background: rgba(63,182,24,.15);  color: var(--success); border: 1px solid rgba(63,182,24,.3); }
.status-err    { background: rgba(255,0,57,.15);    color: var(--danger);  border: 1px solid rgba(255,0,57,.3); }

/* ── Paste textareas ────────────────────────────────────── */
.paste-area         { background: var(--bg) !important; border: 1px solid var(--border) !important; color: var(--text) !important; font-family: var(--mono); font-size: 11px; border-radius: var(--radius); width: 100%; resize: vertical, horizontal; padding: 8px; }
.paste-area:focus   { border-color: var(--primary) !important; outline: none; }
.import-result      { font-family: var(--mono); font-size: 12px; margin-top: 8px; }
.import-ok          { color: var(--success); }
.import-err         { color: var(--danger); }

/* ── Misc helpers ───────────────────────────────────────── */
hr            { border-color: var(--border) !important; }
.hint-text    { color: var(--muted); font-size: 12px; margin-bottom: 10px; }
.section-label { color: var(--muted); font-size: 10px; text-transform: uppercase; letter-spacing: 1px; margin: 12px 0 4px; }
::-webkit-scrollbar        { width: 5px; height: 5px; }
::-webkit-scrollbar-track  { background: var(--bg); }
::-webkit-scrollbar-thumb  { background: var(--border); border-radius: 3px; }
"


# UI ----------------------------------------------------------------------


header <- dashboardHeader(
  title = tags$span(style="font-family:'Source Sans 3','Segoe UI',sans-serif;font-size:16px;font-weight:600;letter-spacing:.2px;","(α)typicalquant Portfolio Optimizer"),
  titleWidth = 325)

sidebar <- dashboardSidebar(width = 255, tags$style(HTML(APP_CSS)),
                            sidebarMenu(id="sidebar",
                                        menuItem("Assets & Assumptions", tabName="assets",       icon=icon("table")),
                                        menuItem("Correlations",         tabName="correlations", icon=icon("th")),
                                        menuItem("Constraints",          tabName="constraints",  icon=icon("sliders-h")),
                                        menuItem("Efficient Frontier",   tabName="frontier",     icon=icon("chart-line")),
                                        menuItem("Portfolio Analysis",   tabName="analysis",     icon=icon("pie-chart")),
                                        menuItem("Portfolio Comparison",   tabName="port_comparison", icon=icon("bar-chart")),
                                        menuItem("Black-Litterman",      tabName="bl",           icon=icon("eye")),
                                        menuItem("Monte Carlo",          tabName="montecarlo",   icon=icon("random")),
                                        menuItem("Comparison & Export",  tabName="comparison",   icon=icon("balance-scale"))))

body <- dashboardBody(useShinyjs(),
                      
                      # Disclaimer Popup
                      
                      tags$div(id="disclaimer_modal",
                               style="display:none; position:fixed; top:0; left:0; width:100%; height:100%;
         background:rgba(0,0,0,0.75); z-index:9999; justify-content:center; align-items:center;",
                               tags$div(
                                 style="background:#1a1f2e; border:1px solid #2d3555; border-radius:4px;
           padding:30px 40px; max-width:560px; margin:auto; margin-top:15vh;",
                                 tags$h4("Important Notice", style="color:#d0d6e8; font-family:'Source Sans 3'; margin-bottom:16px;"),
                                 tags$p("This tool is for educational purposes only and should not be relied upon for making investment decisions.", style="color:#7a87a8; font-size:13px; line-height:1.6;"),
                                 tags$p("The default return, risk, and correlation assumptions are entirely made up and have no basis in reality.", style="color:#7a87a8; font-size:13px; line-height:1.6;"),
                                 tags$p("This tool was generated with the help of AI tools and has not been thouroughly tested.", style="color:#7a87a8; font-size:13px; line-height:1.6;"),
                                 actionButton("dismiss_disclaimer", "I Understand", class="btn-primary",
                                              style="margin-top:16px; width:100%;")
                               )
                      )
                      
                      , tabItems(
  
  ## 1. Assets ---------------------------------------------------------------
  
  tabItem(tabName="assets",
          fluidRow(box(width=12, title="Risk-Free Rate Assumption"
                       # , p(class="hint-text","Input Risk-Free Rate Assumption") 
                       , numericInput("rf_rate","Risk-Free Rate (% p.a.)",value=1.5,min=0,max=15,step=0.1)
          )),
          fluidRow(box(width=12, title="Asset Class Assumptions",
                       p(class="hint-text","Edit returns (%) and volatilities (%) directly in the table. The Class column drives asset-class constraints."),
                       fluidRow(
                         column(3, textInput("new_asset_name","Name",placeholder="e.g., Private Equity")),
                         column(2, numericInput("new_asset_ret","Return (%)",value=5.0,step=0.5)),
                         column(2, numericInput("new_asset_vol","Vol (%)",value=12.0,step=0.5)),
                         column(3, textInput("new_asset_class","Class",placeholder="e.g., Alternatives")),
                         column(2, br(), actionButton("add_asset","Add Asset",class="btn-primary"))),
                       fluidRow(column(12,
                                       actionButton("reset_assets","Reset Defaults",class="btn-default"),
                                       actionButton("remove_asset","Remove Selected",class="btn-danger",style="margin-left:8px;"))),
                       br(), DTOutput("assets_table"))),
          fluidRow(box(width=12, title="Paste Assets from Spreadsheet",
                       p(class="hint-text", "Paste directly from Excel or Google Sheets. Expected columns (with a header row): ",
                         tags$b("Asset, Return (%), Vol (%), Class, Current (%)"), ". Tab- and comma-separated formats are both detected automatically."),
                       tags$textarea(id="paste_assets_raw", class="paste-area", rows="8",
                                     placeholder="Asset\t\tReturn (%)\tVol (%)\tClass\t\tCurrent\nUS Equity\t7.5\t\t16.0\tEquity\t\t60\nUS Bonds\t3.0\t\t5.0\tFixed Income\t40"),
                       br(),
                       fluidRow(
                         column(6, checkboxInput("paste_assets_replace","Replace all existing assets",value=TRUE))),
                       fluidRow(
                         column(6, actionButton("import_assets","Import Assets",class="btn-primary"))),
                       uiOutput("paste_assets_result")
          ))),
  
  ## 2. Correlations ---------------------------------------------------------
  
  tabItem(tabName="correlations",
          fluidRow(box(width=12, title="Correlation Matrix",
                       p(class="hint-text","Pairwise correlations (−100 to +100). Kept symmetric automatically. Diagonal locked at 100."),
                       DTOutput("cor_table"))),
          fluidRow(box(width=12, title="Paste Correlation Matrix from Spreadsheet",
                       p(class="hint-text", "Paste a square matrix from Excel. Include a header row and asset names in the first column — ",
                         "names must match your asset list exactly. Values should be in the −100 to +100 range."),
                       tags$textarea(id="paste_cor_raw", class="paste-area", rows="10",
                                     placeholder="Asset\t\tUS Equity\tUS Bonds\nUS Equity\t100\t\t-10\nUS Bonds\t-10\t\t100"),
                       br(),
                       actionButton("import_cor","Import Correlations",class="btn-primary"),
                       uiOutput("paste_cor_result")
          ))),
  
  ## 3. Constraints ----------------------------------------------------------
  
  tabItem(tabName="constraints",
          fluidRow(
          column(6,
                 fluidRow(
                   # box(width=5, title="Global & Solver Settings",
                   #   sliderInput("w_min_global","Global Min Weight per Asset (%)",min=-30,max=20,value=0,step=1),
                   #   sliderInput("w_max_global","Global Max Weight per Asset (%)",min=5,max=100,value=100,step=5),
                   #   p(class="hint-text","Negative minimums allow short positions. Per-asset overrides below take precedence.")
                   #   # ,
                   #   # hr()
                   #   # ,
                   #   # numericInput("rf_rate","Risk-Free Rate (% p.a.)",value=1.5,min=0,max=15,step=0.1),
                   #   # numericInput("n_frontier","Frontier Points",value=60,min=20,max=200,step=10)
                   #   ),
                   box(width=12, title="Asset Weight Constraints",
                       sliderInput("w_min_global","Global Min Weight per Asset (%)",min=-30,max=20,value=0,step=1),
                       sliderInput("w_max_global","Global Max Weight per Asset (%)",min=5,max=100,value=100,step=5),
                       p(class="hint-text","Negative minimums allow short positions. Per-asset overrides below take precedence."),
                       hr(),
                       p(class="hint-text","Override global min/max for individual assets. Leave blank to inherit the global setting."),
                       uiOutput("per_asset_table_ui")))),
          column(6,
                 fluidRow(
                   box(width=12, title="Asset Class Constraints",
                       p(class="hint-text","Constrain total allocation to an entire asset class (e.g., all Equity between 30–60%)."),
                       uiOutput("class_constraint_ui"), br(),
                       actionButton("add_class_con","Add Class Constraint",class="btn-default"),
                       actionButton("remove_class_con","Remove Last",class="btn-danger",style="margin-left:8px;"))),
                 fluidRow(
                   box(width=12, title="Custom Group Constraints",
                       p(class="hint-text","Combine any assets into a custom group with combined allocation bounds."),
                       uiOutput("group_constraint_ui"), br(),
                       actionButton("add_gc","Add Group Constraint",class="btn-default"),
                       actionButton("remove_gc","Remove Last",class="btn-danger",style="margin-left:8px;"))),
                 fluidRow(box(width=12, title="Solver Settings"
                              , numericInput("n_frontier","Frontier Points",value=60,min=20,max=200,step=10)
                 ))
          ))
          ),
  
  ## 4. Frontier -------------------------------------------------------------
  
  tabItem(tabName="frontier",
          fluidRow(column(12,
                          actionButton("run_optimizer","Run Optimizer",class="btn-primary",
                                       style="font-size:14px;padding:10px 30px;margin-bottom:15px;"),
                          uiOutput("frontier_status"))),
          fluidRow(box(width=12, title="Efficient Frontier", plotlyOutput("frontier_plot",height="500px"))),
          fluidRow(
            box(width=4, title="Global Min-Variance", uiOutput("gmv_stats")),
            box(width=4, title="Max Sharpe Ratio",    uiOutput("msr_stats")),
            box(width=4, title="Selected Portfolio",  uiOutput("sel_stats")))),
  
  ## 5a. Analysis -------------------------------------------------------------
  
  tabItem(tabName="analysis",
          fluidRow(box(width=12, title="Portfolio Selection",
                       fluidRow(
                         column(4, selectInput("selected_objective","Objective Function", choices=c(
                           "Minimum Variance"          = "min_var",
                           "Maximum Sharpe Ratio"      = "max_sharpe",
                           "Maximum Diversification"   = "max_div",
                           "Risk Parity"               = "risk_parity",
                           "Equal Weight"              = "equal_weight",
                           "Min CVaR (5%)"             = "min_cvar",
                           "Maximum Return"            = "max_return",
                           "Maximum Utility"           = "max_utility",
                           "Min Tracking Error"        = "min_te",
                           "Target Return"             = "target_ret",
                           "Target Volatility"         = "target_vol",
                           "Manual Weights"            = "manual_weights",
                           "Current Portfolio"         = "current_port"))),
                         column(4, uiOutput("target_input_ui")),
                         column(4, br(), actionButton("analyse_portfolio","Analyse",class="btn-primary"))))),
          fluidRow(
            valueBoxOutput("vb_return",width=3), valueBoxOutput("vb_vol",width=3),
            valueBoxOutput("vb_sharpe",width=3), valueBoxOutput("vb_cvar",width=3)),
          fluidRow(
            box(width=6, title="Asset Allocation",  plotlyOutput("weights_bar",height="360px")),
            box(width=6, title="Risk Contribution", plotlyOutput("rc_bar",height="360px"))),
          fluidRow(box(width=12, title="Portfolio Weights & Risk", DTOutput("weights_table")))),
  
  ## 5b. Portfolio Comparison ------------------------------------------------
  
  tabItem(tabName="port_comparison",
          fluidRow(box(width=12, title="Compare Portfolios Side by Side",
                       p(class="hint-text","Select up to 3 portfolios to compare directly. Target Return/Volatility slots use the values set in the main Portfolio Selection panel above."),
                       fluidRow(
                         column(4,
                                selectInput("cmp_port1","Portfolio A",choices=c(
                                  "— None —"                  = "none",
                                  "Current Portfolio"         = "current_port",
                                  "Min Variance"              = "min_var",
                                  "Max Sharpe"                = "max_sharpe",
                                  "Max Diversification"       = "max_div",
                                  "Risk Parity"               = "risk_parity",
                                  "Equal Weight"              = "equal_weight",
                                  "Min CVaR (5%)"             = "min_cvar",
                                  "Max Return"                = "max_return",
                                  "Max Utility"               = "max_utility",
                                  "Target Return"             = "target_ret_c",
                                  "Target Volatility"         = "target_vol_c",
                                  "Manual Weights"            = "manual_weights",
                                  "Current Analysis"          = "analysis_result"
                                ),selected="current_port"),
                                uiOutput("cmp_extra1_ui")),
                         column(4,
                                selectInput("cmp_port2","Portfolio B",choices=c(
                                  "— None —"                  = "none",
                                  "Current Portfolio"         = "current_port",
                                  "Min Variance"              = "min_var",
                                  "Max Sharpe"                = "max_sharpe",
                                  "Max Diversification"       = "max_div",
                                  "Risk Parity"               = "risk_parity",
                                  "Equal Weight"              = "equal_weight",
                                  "Min CVaR (5%)"             = "min_cvar",
                                  "Max Return"                = "max_return",
                                  "Max Utility"               = "max_utility",
                                  "Target Return"             = "target_ret_c",
                                  "Target Volatility"         = "target_vol_c",
                                  "Manual Weights"            = "manual_weights",
                                  "Current Analysis"          = "analysis_result"
                                ),selected="max_sharpe"),
                                uiOutput("cmp_extra2_ui")),
                         column(4,
                                selectInput("cmp_port3","Portfolio C",choices=c(
                                  "— None —"                  = "none",
                                  "Current Portfolio"         = "current_port",
                                  "Min Variance"              = "min_var",
                                  "Max Sharpe"                = "max_sharpe",
                                  "Max Diversification"       = "max_div",
                                  "Risk Parity"               = "risk_parity",
                                  "Equal Weight"              = "equal_weight",
                                  "Min CVaR (5%)"             = "min_cvar",
                                  "Max Return"                = "max_return",
                                  "Max Utility"               = "max_utility",
                                  "Target Return"             = "target_ret_c",
                                  "Target Volatility"         = "target_vol_c",
                                  "Manual Weights"            = "manual_weights",
                                  "Current Analysis"          = "analysis_result"
                                ),selected="min_var"),
                                uiOutput("cmp_extra3_ui"))),
                       fluidRow(column(12, br(),
                                       actionButton("compare_portfolios","Compare",class="btn-primary"))),
                       br(),
                       uiOutput("compare_stats_ui"),
                       br(),
                       plotlyOutput("compare_weights_plot", height="360px"),
                       br(),
                       DTOutput("compare_table")
          ))),
  
  ## 6. Black-Litterman ------------------------------------------------------
  
  tabItem(tabName="bl",
          fluidRow(
            box(width=5, title="Market Prior & Settings",
                p(class="hint-text","Black-Litterman blends CAPM-implied equilibrium returns with your subjective views using Bayesian updating. The posterior replaces raw expected returns in the optimizer."),
                p(class="section-label","Market-Cap Weights (%)"),
                uiOutput("bl_mkt_weights_ui"), hr(),
                numericInput("bl_delta","Risk Aversion (delta)",value=2.5,min=0.5,max=10,step=0.5),
                numericInput("bl_tau","Prior Confidence (tau)",value=0.05,min=0.01,max=0.5,step=0.01),
                p(class="hint-text","tau controls how tightly the posterior is pulled toward the equilibrium prior. Smaller tau = stronger prior. Typical: 0.01–0.10."),
                hr(), actionButton("bl_run","Compute BL Posterior",class="btn-primary")),
            box(width=7, title="Investor Views",
                p(class="hint-text","Each view specifies a portfolio (P-vector: +long, −short, 0=excluded), an expected return Q, and a confidence level (0=ignore view, 1=certain)."),
                uiOutput("bl_views_ui"), br(),
                fluidRow(
                  column(6, actionButton("bl_add_view","Add View",class="btn-default")),
                  column(6, actionButton("bl_remove_view","Remove Last",class="btn-danger"))))),
          fluidRow(
            box(width=6, title="Equilibrium vs BL Returns", plotlyOutput("bl_returns_plot",height="380px")),
            box(width=6, title="BL Max-Sharpe Portfolio",   plotlyOutput("bl_weights_plot",height="380px"))),
          fluidRow(box(width=12, title="Return Comparison Table", DTOutput("bl_table")))),
  
  ## 7. Monte Carlo ----------------------------------------------------------
  
  tabItem(tabName="montecarlo",
          fluidRow(
            box(width=4, title="Simulation Settings",
                p(class="hint-text","Simulates annual returns from a multivariate normal distribution over a multi-year horizon. Run the Optimizer or Analysis first."),
                selectInput("mc_portfolio","Portfolio to Simulate",choices=c(
                  "Min Variance"          = "min_var",
                  "Max Sharpe"            = "max_sharpe",
                  "Max Diversification"   = "max_div",
                  "Risk Parity"           = "risk_parity",
                  "Equal Weight"          = "equal_weight",
                  "Current Portfolio"     = "current",
                  "Current Analysis"      = "analysis",
                  "Black-Litterman (MSR)" = "bl_msr")),
                numericInput("mc_horizon","Horizon (years)",value=10,min=1,max=50,step=1),
                numericInput("mc_sims","Number of Simulations",value=1000,min=100,max=10000,step=100),
                numericInput("mc_initial","Initial Portfolio Value ($)",value=1000000,min=1000,step=50000),
                numericInput("mc_contrib","Annual Contribution ($, 0=none)",value=0,min=0,step=5000),
                hr(),
                numericInput("mc_alpha","CVaR Confidence Level (%)",value=5,min=1,max=25,step=1),
                hr(), actionButton("run_mc","Run Simulation",class="btn-primary")),
            box(width=8, title="Wealth Path Simulation", plotlyOutput("mc_paths_plot",height="440px"))),
          fluidRow(
            valueBoxOutput("mc_vb_median",width=3), valueBoxOutput("mc_vb_p10",width=3),
            valueBoxOutput("mc_vb_p90",width=3),   valueBoxOutput("mc_vb_cvar",width=3)),
          fluidRow(
            box(width=6, title="Terminal Wealth Distribution", plotlyOutput("mc_hist_plot",height="340px")),
            box(width=6, title="Annual Return Boxplots",       plotlyOutput("mc_box_plot",height="340px"))),
          fluidRow(
            box(width=12, title="Annualised Return Distribution (CAGR across all paths)",
                p(class="hint-text","Distribution of compound annualised growth rates (CAGRs) across all simulated paths over the full horizon. Vertical lines show the median and 5th/95th percentile outcomes."),
                plotlyOutput("mc_cagr_plot", height="320px"))),
          fluidRow(
            box(width=12, title="Tail Risk Summary",
                p(class="hint-text","Best and worst annualised return and maximum drawdown at the 90%, 95%, and 99% confidence levels across all simulated paths."),
                DTOutput("mc_risk_table")))
          ),
  
  ## 8. Comparison & Export --------------------------------------------------
  
  tabItem(tabName="comparison",
          fluidRow(box(width=12, title="Portfolio Comparison & Export",
                       p(class="hint-text","Compare all systematic portfolios. Run the Optimizer first. BL and Analysis portfolios are included if available."),
                       fluidRow(
                         column(4, actionButton("run_comparison","Generate Comparison",class="btn-primary")),
                         column(8,
                                downloadButton("export_csv",  "Export CSV",      class="btn-default"),
                                downloadButton("export_excel","Export Excel",    class="btn-default",style="margin-left:8px;"),
                                downloadButton("export_bl",   "Export BL CSV",  class="btn-default",style="margin-left:8px;"),
                                downloadButton("export_mc",   "Export MC Paths",class="btn-default",style="margin-left:8px;"))))),
          fluidRow(box(width=12, title="Summary Statistics", DTOutput("comparison_stats_table"))),
          fluidRow(box(width=12, title="Weights Heatmap", plotlyOutput("comparison_heatmap",height="460px"))))
  
)) # end tabItems, dashboardBody

ui <- dashboardPage(header, sidebar, body, skin="blue")



# Server ------------------------------------------------------------------


# Auto-detects tab or comma separator; returns list(ok, data/msg)
parse_pasted_table <- function(raw_text) {
  lines <- strsplit(trimws(raw_text), "\n")[[1]]
  lines <- lines[nchar(trimws(lines)) > 0]
  if (length(lines) < 2)
    return(list(ok=FALSE, msg="Need at least a header row and one data row."))
  sep    <- if (grepl("\t", lines[1])) "\t" else ","
  parsed <- lapply(lines, function(l) strsplit(l, sep, fixed=TRUE)[[1]])
  ncols  <- length(parsed[[1]])
  if (!all(sapply(parsed, length) == ncols))
    return(list(ok=FALSE, msg="Inconsistent column count — check for stray separators."))
  df           <- as.data.frame(do.call(rbind, parsed[-1]), stringsAsFactors=FALSE)
  colnames(df) <- trimws(parsed[[1]])
  list(ok=TRUE, data=df)
}

# ──────────────────────────────────────────────────────────────
server <- function(input, output, session) {
  
  shinyjs::runjs("document.title = '(α)typicalquant Portfolio Optimizer';")
  
  # Show disclaimer on launch
  shinyjs::runjs("
  var el = document.getElementById('disclaimer_modal');
  el.style.display = 'flex';
  document.body.style.overflow = 'hidden';
")

  observeEvent(input$dismiss_disclaimer, {
    shinyjs::runjs("
    document.getElementById('disclaimer_modal').style.display = 'none';
    document.body.style.overflow = '';
  ")
  })
  
  rv <- reactiveValues(
    assets=default_assets, cor_matrix=default_cor,
    frontier=NULL, gmv_port=NULL, msr_port=NULL, analysis_w=NULL,
    gc_count=0, class_con_count=0, bl_view_count=0,
    bl_mu=NULL, bl_Sigma=NULL, bl_weights=NULL,
    mc_result=NULL, comp_data=NULL)
  
  ## Shared Reactives --------------------------------------------------------
  
  get_mu_sigma <- reactive({
    mu <- rv$assets$Return_pct/100; vols <- rv$assets$Vol_pct
    list(mu=mu, Sigma=build_cov_matrix(vols, rv$cor_matrix), vols=vols)
  })
  
  get_weight_bounds <- reactive({
    n <- nrow(rv$assets)
    wl <- rep(input$w_min_global/100, n)
    wu <- rep(input$w_max_global/100, n)
    for (i in seq_len(n)) {
      vm <- input[[paste0("pa_min_",i)]]; vx <- input[[paste0("pa_max_",i)]]
      if (!is.null(vm) && !is.na(vm) && nchar(trimws(as.character(vm)))>0) wl[i] <- as.numeric(vm)/100
      if (!is.null(vx) && !is.na(vx) && nchar(trimws(as.character(vx)))>0) wu[i] <- as.numeric(vx)/100
    }
    list(w_min_vec=wl, w_max_vec=wu)
  })
  
  get_all_constraints <- reactive({
    assets <- rv$assets; gc_list <- list()
    if (rv$class_con_count > 0) {
      for (i in seq_len(rv$class_con_count)) {
        cls <- input[[paste0("cc_class_",i)]]; mn <- input[[paste0("cc_min_",i)]]; mx <- input[[paste0("cc_max_",i)]]
        if (is.null(cls)||cls=="") next
        idx <- which(assets$Class==cls); if (length(idx)==0) next
        gc_list <- c(gc_list, list(list(indices=idx, min=mn%||%0, max=mx%||%100)))
      }
    }
    if (rv$gc_count > 0) {
      for (i in seq_len(rv$gc_count)) {
        sel <- input[[paste0("gc_assets_",i)]]; mn <- input[[paste0("gc_min_",i)]]; mx <- input[[paste0("gc_max_",i)]]
        if (is.null(sel)||length(sel)==0) next
        idx <- which(assets$Asset %in% sel)
        gc_list <- c(gc_list, list(list(indices=idx, min=mn%||%0, max=mx%||%100)))
      }
    }
    if (length(gc_list)==0) NULL else gc_list
  })
  
  # Returns normalised current portfolio weight vector, or NULL if no weights entered
  get_current_portfolio <- reactive({
    w_raw <- rv$assets$Current_pct
    if (is.null(w_raw) || sum(w_raw, na.rm=TRUE) < 0.001) return(NULL)
    w_raw[is.na(w_raw)] <- 0
    w_raw / sum(w_raw)
  })
  
  plot_dark <- function(p, xlab="", ylab="", legend=FALSE) {
    p %>% layout(
      paper_bgcolor="#1a1f2e", plot_bgcolor="#0f1117",
      xaxis=list(title=list(text=xlab,font=list(color="#7a87a8",family="Source Sans 3",size=12)),
                 tickfont=list(color="#7a87a8",family="IBM Plex Mono",size=11),gridcolor="#2d3555",zerolinecolor="#2d3555"),
      yaxis=list(title=list(text=ylab,font=list(color="#7a87a8",family="Source Sans 3",size=12)),
                 tickfont=list(color="#7a87a8",family="IBM Plex Mono",size=11),gridcolor="#2d3555",zerolinecolor="#2d3555"),
      legend=list(font=list(color="#7a87a8",family="Source Sans 3",size=12),bgcolor="#1a1f2e",bordercolor="#2d3555"),
      showlegend=legend, margin=list(l=60,r=15,t=15,b=70))
  }
  
  ## 1. Assets ---------------------------------------------------------------
  
  
  output$assets_table <- renderDT({
    datatable(rv$assets, editable=list(target="cell"), selection="single",
              rownames=FALSE, class="compact stripe",
              colnames=c("Asset","Return (%)","Volatility (%)","Class","Current (%)"),
              options=list(dom="t",paging=FALSE,scrollX=TRUE,
                           columnDefs=list(list(className="dt-center",targets=4))))
  }, server=FALSE)
  
  observeEvent(input$assets_table_cell_edit, {
    info <- input$assets_table_cell_edit; df <- rv$assets; ci <- info$col+1
    df[info$row,ci] <- if (ci %in% c(2,3,5)) as.numeric(info$value) else as.character(info$value)
    rv$assets <- df
    n_new <- nrow(df); n_old <- nrow(rv$cor_matrix)
    if (n_new!=n_old) {
      nc <- matrix(30,n_new,n_new); diag(nc) <- 100; mn <- min(n_old,n_new)
      nc[1:mn,1:mn] <- rv$cor_matrix[1:mn,1:mn]; rv$cor_matrix <- nc
    }
  })
  
  observeEvent(input$add_asset, {
    nm <- trimws(input$new_asset_name)
    if (nchar(nm)==0) { showNotification("Enter an asset name.",type="warning"); return() }
    cls <- trimws(input$new_asset_class); if (nchar(cls)==0) cls <- "Other"
    rv$assets <- rbind(rv$assets,
                       data.frame(Asset=nm,Return_pct=input$new_asset_ret,Vol_pct=input$new_asset_vol,
                                  Class=cls,Current_pct=0,stringsAsFactors=FALSE))
    n <- nrow(rv$assets)
    rv$cor_matrix <- rbind(cbind(rv$cor_matrix,30),c(rep(30,n-1),100))
    updateTextInput(session,"new_asset_name",value="")
    updateTextInput(session,"new_asset_class",value="")
  })
  
  observeEvent(input$remove_asset, {
    sel <- input$assets_table_rows_selected
    if (is.null(sel)) { showNotification("Select a row first.",type="warning"); return() }
    if (nrow(rv$assets)<=2) { showNotification("Need at least 2 assets.",type="warning"); return() }
    rv$assets <- rv$assets[-sel,]; rv$cor_matrix <- rv$cor_matrix[-sel,-sel]
  })
  
  observeEvent(input$reset_assets, {
    rv$assets <- default_assets; rv$cor_matrix <- default_cor
    showNotification("Defaults restored.",type="message")
  })
  
  # Pre-populate Min Tracking Error benchmark inputs from the Current Portfolio
  observeEvent(input$bench_from_current, {
    n     <- nrow(rv$assets)
    w_cur <- rv$assets$Current_pct
    if (sum(w_cur, na.rm=TRUE) < 0.001) return()
    for (i in seq_len(n))
      updateNumericInput(session, paste0("bench_w_",i), value=round(w_cur[i],1))
  })
  
  observeEvent(input$import_assets, {
    raw <- input$paste_assets_raw
    if (is.null(raw) || nchar(trimws(raw)) == 0) {
      output$paste_assets_result <- renderUI(p(class="import-result import-err","Nothing pasted yet.")); return()
    }
    res <- parse_pasted_table(raw)
    if (!res$ok) {
      output$paste_assets_result <- renderUI(p(class="import-result import-err", paste("Parse error:", res$msg))); return()
    }
    df   <- res$data
    cols <- tolower(colnames(df))
    find_col <- function(patterns) {
      idx <- which(sapply(cols, function(c) any(sapply(patterns, function(p) grepl(p, c, fixed=FALSE)))))
      if (length(idx) == 0) NA_integer_ else idx[1]
    }
    ci_name  <- find_col(c("asset","name","ticker","symbol"))
    ci_ret   <- find_col(c("return","ret","mu","expected"))
    ci_vol   <- find_col(c("vol","risk","sigma","std","sd"))
    ci_class <- find_col(c("class","type","category","group"))
    if (is.na(ci_name) || is.na(ci_ret) || is.na(ci_vol)) {
      output$paste_assets_result <- renderUI(p(class="import-result import-err",
                                               "Could not identify required columns. Need: Asset name, Return (%), Volatility (%).")); return()
    }
    tryCatch({
      ci_cur <- find_col(c("current","existing","actual","portfolio","weight","alloc"))
      new_df <- data.frame(
        Asset       = trimws(df[[ci_name]]),
        Return_pct  = as.numeric(df[[ci_ret]]),
        Vol_pct     = as.numeric(df[[ci_vol]]),
        Class       = if (!is.na(ci_class)) trimws(df[[ci_class]]) else rep("Other", nrow(df)),
        Current_pct = if (!is.na(ci_cur))  as.numeric(df[[ci_cur]]) else rep(0, nrow(df)),
        stringsAsFactors = FALSE
      )
      if (any(is.na(new_df$Return_pct)) || any(is.na(new_df$Vol_pct)))
        stop("Non-numeric values found in Return or Volatility columns.")
      if (nrow(new_df) < 2) stop("Need at least 2 assets.")
      n_new <- nrow(new_df)
      if (isTRUE(input$paste_assets_replace)) {
        rv$assets     <- new_df
        new_cor       <- matrix(30, n_new, n_new); diag(new_cor) <- 100
        rv$cor_matrix <- new_cor
      } else {
        rv$assets <- rbind(rv$assets, new_df)
        n_tot <- nrow(rv$assets)
        new_cor <- matrix(30, n_tot, n_tot); diag(new_cor) <- 100
        n_old <- nrow(rv$cor_matrix)
        new_cor[1:n_old, 1:n_old] <- rv$cor_matrix
        rv$cor_matrix <- new_cor
      }
      output$paste_assets_result <- renderUI(
        p(class="import-result import-ok", paste0("✓ ", n_new, " assets imported successfully.")))
    }, error = function(e) {
      output$paste_assets_result <- renderUI(
        p(class="import-result import-err", paste("Import error:", e$message)))
    })
  })
  
  observeEvent(input$import_cor, {
    raw <- input$paste_cor_raw
    if (is.null(raw) || nchar(trimws(raw)) == 0) {
      output$paste_cor_result <- renderUI(p(class="import-result import-err","Nothing pasted yet.")); return()
    }
    res <- parse_pasted_table(raw)
    if (!res$ok) {
      output$paste_cor_result <- renderUI(p(class="import-result import-err", paste("Parse error:", res$msg))); return()
    }
    df <- res$data; n <- nrow(rv$assets)
    tryCatch({
      row_labels <- trimws(df[[1]])
      mat_raw    <- df[, -1, drop=FALSE]
      if (ncol(mat_raw) != nrow(mat_raw))
        stop("Matrix is not square after removing the label column.")
      mat <- apply(mat_raw, 2, as.numeric)
      if (any(is.na(mat))) stop("Non-numeric values found in the matrix.")
      if (nrow(mat) != n) stop(paste0(
        "Matrix is ", nrow(mat), "×", ncol(mat), " but app has ", n,
        " assets. Import or add assets first if needed."))
      if (any(abs(mat) > 100)) stop("Values outside −100 to +100 detected.")
      # Reorder to match current asset order if names are present and matchable
      asset_names <- rv$assets$Asset
      col_labels  <- trimws(colnames(mat_raw))
      row_match   <- match(asset_names, row_labels)
      col_match   <- match(asset_names, col_labels)
      if (!anyNA(row_match) && !anyNA(col_match)) mat <- mat[row_match, col_match]
      # Enforce symmetry and lock diagonal
      mat <- (mat + t(mat)) / 2; diag(mat) <- 100
      rv$cor_matrix <- mat
      output$paste_cor_result <- renderUI(
        p(class="import-result import-ok", paste0("✓ ", n, "×", n, " correlation matrix imported successfully.")))
    }, error = function(e) {
      output$paste_cor_result <- renderUI(
        p(class="import-result import-err", paste("Import error:", e$message)))
    })
  })
  
  ## 2. Correlations ---------------------------------------------------------
  
  
  output$cor_table <- renderDT({
    df <- as.data.frame(rv$cor_matrix); colnames(df) <- rv$assets$Asset
    datatable(df, editable=list(target="cell"), rownames=rv$assets$Asset,
              class="compact stripe", options=list(dom="t",paging=FALSE,scrollX=TRUE))
  }, server=FALSE)
  
  observeEvent(input$cor_table_cell_edit, {
    info <- input$cor_table_cell_edit; r <- info$row; c <- info$col
    val <- max(-100,min(100,as.numeric(info$value))); if (r==c) val <- 100
    rv$cor_matrix[r,c] <- val; rv$cor_matrix[c,r] <- val
  })
  
  ## 3. Constraints ----------------------------------------------------------
  
  output$per_asset_table_ui <- renderUI({
    n <- nrow(rv$assets)
    tagList(tags$table(style="width:100%;border-collapse:collapse;",
                       tags$thead(tags$tr(
                         tags$th(style="color:#d0d6e8;font-family:'Source Sans 3',sans-serif;font-size:11px;text-transform:uppercase;letter-spacing:.7px;padding:6px 8px;background:#222840;width:45%;","Asset"),
                         tags$th(style="color:#d0d6e8;font-family:'Source Sans 3',sans-serif;font-size:11px;text-transform:uppercase;letter-spacing:.7px;padding:6px 8px;background:#222840;","Min (%)"),
                         tags$th(style="color:#d0d6e8;font-family:'Source Sans 3',sans-serif;font-size:11px;text-transform:uppercase;letter-spacing:.7px;padding:6px 8px;background:#222840;","Max (%)"))),
                       tags$tbody(lapply(seq_len(n), function(i) {
                         tags$tr(style=if(i%%2==0)"background:#222840;" else "background:#1a1f2e;",
                                 tags$td(style="color:#d0d6e8;font-family:'IBM Plex Mono';font-size:12px;padding:4px 8px;",rv$assets$Asset[i]),
                                 tags$td(style="padding:3px 5px;",numericInput(paste0("pa_min_",i),NULL,value=NA,min=-100,max=100,step=1,width="100%")),
                                 tags$td(style="padding:3px 5px;",numericInput(paste0("pa_max_",i),NULL,value=NA,min=0,max=100,step=1,width="100%")))
                       }))))
  })
  
  output$class_constraint_ui <- renderUI({
    if (rv$class_con_count==0) return(p(class="hint-text","No class constraints defined."))
    classes <- unique(rv$assets$Class)
    lapply(seq_len(rv$class_con_count), function(i) {
      wellPanel(class="well-sm", fluidRow(
        column(4, selectInput(paste0("cc_class_",i),paste("Constraint",i,"Class"),
                              choices=classes,selected=classes[min(i,length(classes))])),
        column(4, numericInput(paste0("cc_min_",i),"Min (%)",value=0,min=0,max=100)),
        column(4, numericInput(paste0("cc_max_",i),"Max (%)",value=60,min=0,max=100))))
    })
  })
  observeEvent(input$add_class_con,    { rv$class_con_count <- rv$class_con_count+1 })
  observeEvent(input$remove_class_con, { if (rv$class_con_count>0) rv$class_con_count <- rv$class_con_count-1 })
  
  output$group_constraint_ui <- renderUI({
    if (rv$gc_count==0) return(p(class="hint-text","No custom group constraints defined."))
    lapply(seq_len(rv$gc_count), function(i) {
      wellPanel(class="well-sm",
                tags$b(style="color:var(--text);font-size:13px;font-weight:600;",paste("Group",i)),
                checkboxGroupInput(paste0("gc_assets_",i),"Assets:",choices=rv$assets$Asset,inline=TRUE),
                fluidRow(
                  column(6,numericInput(paste0("gc_min_",i),"Min (%)",value=0,min=0,max=100)),
                  column(6,numericInput(paste0("gc_max_",i),"Max (%)",value=60,min=0,max=100))))
    })
  })
  observeEvent(input$add_gc,    { rv$gc_count <- rv$gc_count+1 })
  observeEvent(input$remove_gc, { if (rv$gc_count>0) rv$gc_count <- rv$gc_count-1 })
  
  ### Optimizer Core ----------------------------------------------------------
  
  observeEvent(input$run_optimizer, {
    req(nrow(rv$assets)>=2)
    ms <- get_mu_sigma(); wb <- get_weight_bounds(); gc <- get_all_constraints()
    withProgress(message="Optimising...",value=0,{
      setProgress(0.1,detail="Generating frontier")
      front <- generate_frontier(ms$mu,ms$Sigma,n_points=input$n_frontier,
                                 w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc)
      setProgress(0.9,detail="Key portfolios")
      if (is.null(front)) { showNotification("Frontier failed — constraints may be infeasible.",type="error"); return() }
      rv$frontier <- front; rv$gmv_port <- front[[1]]
      rv$msr_port <- max_sharpe_portfolio(front,rf=input$rf_rate/100)
      rv$analysis_w <- NULL; setProgress(1)
    })
    showNotification(paste0("Frontier: ",length(rv$frontier)," feasible portfolios."),type="message")
  })
  
  ## 4. Frontier -------------------------------------------------------------
  
  output$frontier_status <- renderUI({
    if (is.null(rv$frontier))
      div(style="padding:8px 0;margin-bottom:10px;",span(class="status-badge status-err","Not computed — click Run Optimizer"))
    else
      div(style="padding:8px 0;margin-bottom:10px;",span(class="status-badge status-ok",paste0("✓ ",length(rv$frontier)," feasible portfolios on frontier")))
  })
  
  output$frontier_plot <- renderPlotly({
    req(rv$frontier); ms <- get_mu_sigma()
    vols <- sapply(rv$frontier,`[[`,"Vol"); rets <- sapply(rv$frontier,`[[`,"Return")*100
    gmv_s <- port_stats(rv$gmv_port$Weights,ms$mu,ms$Sigma,input$rf_rate/100)
    msr_s <- port_stats(rv$msr_port$Weights,ms$mu,ms$Sigma,input$rf_rate/100)
    p <- plot_ly() %>%
      add_trace(x=vols,y=rets,type="scatter",mode="lines",
                line=list(color="#2780e3",width=2.5),name="Efficient Frontier",
                hovertemplate="Vol: %{x:.2f}%<br>Return: %{y:.2f}%<extra></extra>") %>%
      add_trace(x=rv$assets$Vol_pct,y=rv$assets$Return_pct,type="scatter",mode="markers+text",
                marker=list(color="#3fb618",size=9,symbol="diamond",line=list(color="#d0d6e8",width=1)),
                text=rv$assets$Asset,textposition="top right",
                textfont=list(family="Source Sans 3",size=10,color="#7a87a8"),name="Assets",
                hovertemplate="%{text}<br>Vol: %{x:.1f}%<br>Return: %{y:.1f}%<extra></extra>") %>%
      add_trace(x=gmv_s$Vol,y=gmv_s$Return*100,type="scatter",mode="markers+text",
                marker=list(color="#3fb618",size=13,symbol="circle",line=list(color="#d0d6e8",width=1.5)),
                text="GMV",textposition="top left",textfont=list(family="IBM Plex Mono",size=11,color="#3fb618"),
                name="Global Min-Var",hovertemplate="GMV<br>Vol: %{x:.2f}%<br>Return: %{y:.2f}%<extra></extra>") %>%
      add_trace(x=msr_s$Vol,y=msr_s$Return*100,type="scatter",mode="markers+text",
                marker=list(color="#ff0039",size=13,symbol="star",line=list(color="#d0d6e8",width=1.5)),
                text="MSR",textposition="top right",textfont=list(family="IBM Plex Mono",size=11,color="#ff0039"),
                name="Max Sharpe",hovertemplate="MSR<br>Vol: %{x:.2f}%<br>Return: %{y:.2f}%<extra></extra>")
    if (!is.null(rv$analysis_w)) {
      as_s <- port_stats(rv$analysis_w,ms$mu,ms$Sigma,input$rf_rate/100)
      p <- p %>% add_trace(x=as_s$Vol,y=as_s$Return*100,type="scatter",mode="markers+text",
                           marker=list(color="#9954bb",size=12,symbol="triangle-up",line=list(color="#d0d6e8",width=1.5)),
                           text="Selected",textposition="bottom right",textfont=list(family="IBM Plex Mono",size=11,color="#9954bb"),
                           name="Selected Portfolio")
    }
    w_cur <- get_current_portfolio()
    if (!is.null(w_cur)) {
      cur_s <- port_stats(w_cur,ms$mu,ms$Sigma,input$rf_rate/100)
      p <- p %>% add_trace(x=cur_s$Vol,y=cur_s$Return*100,type="scatter",mode="markers+text",
                           marker=list(color="#f5a623",size=12,symbol="square",line=list(color="#d0d6e8",width=1.5)),
                           text="Current",textposition="top left",textfont=list(family="IBM Plex Mono",size=11,color="#f5a623"),
                           name="Current Portfolio")
    }
    plot_dark(p,"Annualised Volatility (%)","Expected Return (%)",legend=TRUE)
  })
  
  stats_html <- function(pstats) {
    rows <- list(
      list("Return",     paste0(round(pstats$Return*100,2),"%"), "#d0d6e8"),
      list("Volatility", paste0(round(pstats$Vol,2),"%"),        "#d0d6e8"),
      list("Sharpe",     round(pstats$Sharpe%||%NA,3),           "#2780e3"))
    tags$table(style="width:100%;",
               do.call(tagList, lapply(rows, function(r)
                 tags$tr(
                   tags$td(style="color:#7a87a8;font-family:'Source Sans 3',sans-serif;font-size:12px;padding:3px 0;",r[[1]]),
                   tags$td(style=paste0("color:",r[[3]],";font-family:'IBM Plex Mono',monospace;font-size:14px;text-align:right;"),r[[2]])))))
  }
  output$gmv_stats <- renderUI({ req(rv$gmv_port); ms<-get_mu_sigma(); stats_html(port_stats(rv$gmv_port$Weights,ms$mu,ms$Sigma,input$rf_rate/100)) })
  output$msr_stats <- renderUI({ req(rv$msr_port); ms<-get_mu_sigma(); stats_html(port_stats(rv$msr_port$Weights,ms$mu,ms$Sigma,input$rf_rate/100)) })
  output$sel_stats <- renderUI({ req(rv$analysis_w); ms<-get_mu_sigma(); stats_html(port_stats(rv$analysis_w,ms$mu,ms$Sigma,input$rf_rate/100)) })
  
  ## 5. Analysis -------------------------------------------------------------
  
  output$target_input_ui <- renderUI({
    switch(input$selected_objective,
           "target_ret"  = numericInput("target_return_val","Target Return (% p.a.)",value=5.0,step=0.5),
           "target_vol"  = numericInput("target_vol_val","Target Volatility (%)",value=10.0,step=0.5),
           "max_utility" = numericInput("lambda_val","Risk Aversion (lambda)",value=3.0,min=0.5,max=20,step=0.5),
           "min_te" = {
             n <- nrow(rv$assets); eq <- round(100/n,1)
             w_cur <- rv$assets$Current_pct
             has_cur <- !is.null(w_cur) && sum(w_cur, na.rm=TRUE) > 0.001
             tagList(
               p(class="hint-text","Benchmark weights (%). Click the button to seed from the Current Portfolio."),
               if (has_cur) actionButton("bench_from_current","Use Current Portfolio as Benchmark",
                                         class="btn-default",style="margin-bottom:8px;"),
               lapply(seq_len(n), function(i)
                 numericInput(paste0("bench_w_",i),rv$assets$Asset[i],
                              value=if(has_cur) round(w_cur[i],1) else eq,
                              min=0,max=100,step=0.5,width="100%")))
           },
           "manual_weights" = {
             n  <- nrow(rv$assets); eq <- round(100/n, 1)
             tagList(
               p(class="hint-text","Enter weights (%). Negative values are allowed for short positions. Weights are normalised to sum to 100% automatically."),
               lapply(seq_len(n), function(i)
                 numericInput(paste0("manual_w_",i), rv$assets$Asset[i], value=eq, min=-100, max=100, step=0.5, width="100%")),
               textOutput("manual_w_sum_display")
             )
           },
           "current_port" = {
             w_cur <- rv$assets$Current_pct
             if (is.null(w_cur) || sum(w_cur, na.rm=TRUE) < 0.001)
               p(class="hint-text", "No current portfolio weights entered yet. Add them in the Assets tab.")
             else {
               w_norm <- w_cur / sum(w_cur)
               tagList(
                 p(class="hint-text", "Displaying the Current Portfolio entered in the Assets tab."),
                 tags$table(style="width:100%;",
                            lapply(seq_len(nrow(rv$assets)), function(i)
                              tags$tr(
                                tags$td(style="color:#7a87a8;font-size:12px;padding:2px 0;",rv$assets$Asset[i]),
                                tags$td(style="color:#d0d6e8;font-family:'IBM Plex Mono';font-size:12px;text-align:right;",
                                        paste0(round(w_norm[i]*100,1),"%"))
                              )
                            )
                 )
               )
             }
           },
           NULL)
  })
  
  output$manual_w_sum_display <- renderText({
    req(input$selected_objective == "manual_weights")
    n   <- nrow(rv$assets)
    raw <- sapply(seq_len(n), function(i) input[[paste0("manual_w_",i)]] %||% 0)
    s   <- round(sum(raw), 1)
    paste0("Sum: ", s, "% (will be normalised to 100%)")
  })
  
  observeEvent(input$analyse_portfolio, {
    ms <- get_mu_sigma(); gc <- get_all_constraints(); wb <- get_weight_bounds()
    obj <- input$selected_objective
    w <- switch(obj,
                "min_var"      = solve_mv_portfolio(ms$mu,ms$Sigma,w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc),
                "max_sharpe"   = { req(rv$msr_port); rv$msr_port$Weights },
                "max_div"      = max_div_portfolio(ms$Sigma,ms$vols,wb$w_min_vec,wb$w_max_vec,gc),
                "risk_parity"  = risk_parity_portfolio(ms$Sigma),
                "equal_weight" = equal_weight_portfolio(nrow(rv$assets)),
                "min_cvar"     = withProgress(message="Computing Min-CVaR...",value=0.5,{
                  min_cvar_portfolio(ms$mu,ms$Sigma,alpha=0.05,w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc)}),
                "max_return"   = max_return_portfolio(ms$mu,wb$w_min_vec,wb$w_max_vec),
                "max_utility"  = max_utility_portfolio(ms$mu,ms$Sigma,lambda=input$lambda_val%||%3,w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc),
                "min_te"       = {
                  n <- nrow(rv$assets)
                  wb2 <- sapply(seq_len(n),function(i) input[[paste0("bench_w_",i)]]%||%(100/n))
                  min_te_portfolio(ms$mu,ms$Sigma,wb2/sum(wb2),wb$w_min_vec,wb$w_max_vec,gc)},
                "target_ret"   = { req(input$target_return_val)
                  solve_mv_portfolio(ms$mu,ms$Sigma,target_return=input$target_return_val/100,
                                     w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc) },
                "target_vol"   = { req(input$target_vol_val,rv$frontier)
                  vf <- sapply(rv$frontier,`[[`,"Vol")
                  rv$frontier[[which.min(abs(vf-input$target_vol_val))]]$Weights },
                "manual_weights" = {
                  n   <- nrow(rv$assets)
                  raw <- sapply(seq_len(n), function(i) input[[paste0("manual_w_",i)]] %||% 0)
                  if (sum(abs(raw)) < 1e-6) {
                    showNotification("All weights are zero — enter at least one non-zero weight.", type="warning")
                    return()
                  }
                  raw / sum(raw)
                },
                "current_port" = {
                  w_cur <- get_current_portfolio()
                  if (is.null(w_cur)) {
                    showNotification("No current portfolio weights set — add them in the Assets tab.", type="warning")
                    return()
                  }
                  w_cur
                }
    )
    if (is.null(w)) { showNotification("Optimisation failed — check constraints.",type="error"); return() }
    rv$analysis_w <- w
  })
  
  make_port_df <- reactive({
    req(rv$analysis_w); ms <- get_mu_sigma(); w <- rv$analysis_w
    rc <- risk_contributions(w,ms$Sigma)
    data.frame(Asset=rv$assets$Asset,Class=rv$assets$Class,
               Weight=round(w*100,2),RC_pct=round(rc$RC_pct,2),stringsAsFactors=FALSE)
  })
  
  output$vb_return <- renderValueBox({
    req(rv$analysis_w); ms<-get_mu_sigma(); s<-port_stats(rv$analysis_w,ms$mu,ms$Sigma,input$rf_rate/100)
    valueBox(paste0(round(s$Return*100,2),"%"),"Expected Return",color="yellow",icon=icon("arrow-trend-up"))
  })
  output$vb_vol <- renderValueBox({
    req(rv$analysis_w); ms<-get_mu_sigma(); s<-port_stats(rv$analysis_w,ms$mu,ms$Sigma,input$rf_rate/100)
    valueBox(paste0(round(s$Vol,2),"%"),"Volatility",color="blue",icon=icon("wave-square"))
  })
  output$vb_sharpe <- renderValueBox({
    req(rv$analysis_w); ms<-get_mu_sigma(); s<-port_stats(rv$analysis_w,ms$mu,ms$Sigma,input$rf_rate/100)
    valueBox(round(s$Sharpe%||%NA,3),"Sharpe Ratio",color="green",icon=icon("star"))
  })
  output$vb_cvar <- renderValueBox({
    req(rv$analysis_w); ms<-get_mu_sigma(); set.seed(42)
    sims <- as.numeric(MASS::mvrnorm(5000,mu=ms$mu,Sigma=ms$Sigma) %*% rv$analysis_w)
    q <- quantile(sims,0.05); cvar <- -mean(sims[sims<=q])*100
    valueBox(paste0(round(cvar,2),"%"),"CVaR (5%)",color="red",icon=icon("exclamation-triangle"))
  })
  
  output$weights_bar <- renderPlotly({
    df <- make_port_df(); df <- df[df$Weight>0.01,]
    p <- plot_ly(df,x=~Asset,y=~Weight,type="bar",
                 marker=list(color=PAL[seq_len(nrow(df))],line=list(color="#0f1117",width=1)),
                 hovertemplate="%{x}: %{y:.2f}%<extra></extra>")
    plot_dark(p,ylab="Weight (%)")
  })
  output$rc_bar <- renderPlotly({
    df <- make_port_df(); df <- df[df$Weight>0.01,]
    p <- plot_ly(df,x=~Asset,y=~RC_pct,type="bar",
                 marker=list(color=PAL[seq_len(nrow(df))],opacity=0.75,line=list(color="#0f1117",width=1)),
                 hovertemplate="%{x}: %{y:.2f}%<extra></extra>")
    plot_dark(p,ylab="Risk Contribution (%)")
  })
  output$weights_table <- renderDT({
    df <- make_port_df(); colnames(df) <- c("Asset","Class","Weight (%)","Risk Contribution (%)")
    datatable(df,rownames=FALSE,class="compact stripe",options=list(dom="t",paging=FALSE,scrollX=TRUE)) %>%
      formatStyle("Weight (%)",background=styleColorBar(c(0,100),"rgba(39,128,227,0.2)"),
                  backgroundSize="95% 50%",backgroundRepeat="no-repeat",backgroundPosition="center")
  })
  
  # ── Portfolio comparison helpers ─────────────────────────────
  
  # Extra UI for target-return / target-vol comparison slots
  make_cmp_extra_ui <- function(slot_id, sel) {
    if (sel == "target_ret_c")
      numericInput(paste0("cmp_tret_",slot_id), "Target Return (% p.a.)", value=5.0, step=0.5)
    else if (sel == "target_vol_c")
      numericInput(paste0("cmp_tvol_",slot_id), "Target Volatility (%)", value=10.0, step=0.5)
    else NULL
  }
  output$cmp_extra1_ui <- renderUI({ make_cmp_extra_ui(1, input$cmp_port1 %||% "none") })
  output$cmp_extra2_ui <- renderUI({ make_cmp_extra_ui(2, input$cmp_port2 %||% "none") })
  output$cmp_extra3_ui <- renderUI({ make_cmp_extra_ui(3, input$cmp_port3 %||% "none") })
  
  # Resolve a slot selection to a named weight vector (NULL if none/unavailable)
  resolve_cmp_slot <- function(sel, slot_id) {
    ms  <- get_mu_sigma()
    wb  <- get_weight_bounds()
    gc  <- get_all_constraints()
    switch(sel,
           "none"            = NULL,
           "current_port"    = get_current_portfolio(),
           "analysis_result" = rv$analysis_w,
           "min_var"         = solve_mv_portfolio(ms$mu,ms$Sigma,w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc),
           "max_sharpe"      = { req(rv$msr_port); rv$msr_port$Weights },
           "max_div"         = max_div_portfolio(ms$Sigma,ms$vols,wb$w_min_vec,wb$w_max_vec,gc),
           "risk_parity"     = risk_parity_portfolio(ms$Sigma),
           "equal_weight"    = equal_weight_portfolio(nrow(rv$assets)),
           "min_cvar"        = withProgress(message="Computing Min-CVaR...",value=0.5,{
             min_cvar_portfolio(ms$mu,ms$Sigma,alpha=0.05,
                                w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc)}),
           "max_return"      = max_return_portfolio(ms$mu,wb$w_min_vec,wb$w_max_vec),
           "max_utility"     = max_utility_portfolio(ms$mu,ms$Sigma,lambda=3,
                                                     w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc),
           "target_ret_c"    = {
             tv <- input[[paste0("cmp_tret_",slot_id)]] %||% input$target_return_val %||% 5.0
             solve_mv_portfolio(ms$mu,ms$Sigma,target_return=tv/100,
                                w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc) },
           "target_vol_c"    = {
             req(rv$frontier)
             tv <- input[[paste0("cmp_tvol_",slot_id)]] %||% input$target_vol_val %||% 10.0
             vf <- sapply(rv$frontier,`[[`,"Vol")
             rv$frontier[[which.min(abs(vf - tv))]]$Weights },
           "manual_weights"  = {
             n   <- nrow(rv$assets)
             raw <- sapply(seq_len(n), function(i) input[[paste0("manual_w_",i)]] %||% 0)
             if (sum(abs(raw)) < 1e-6) NULL else raw / sum(raw) },
           NULL
    )
  }
  
  # Labels for each slot selection
  cmp_label <- function(sel) {
    switch(sel,
           "none"="—","current_port"="Current Portfolio","analysis_result"="Current Analysis",
           "min_var"="Min Variance","max_sharpe"="Max Sharpe","max_div"="Max Diversification",
           "risk_parity"="Risk Parity","equal_weight"="Equal Weight","min_cvar"="Min CVaR",
           "max_return"="Max Return","max_utility"="Max Utility",
           "target_ret_c"="Target Return","target_vol_c"="Target Volatility",
           "manual_weights"="Manual Weights", sel)
  }
  
  # Slot accent colours: orange / blue / green
  CMP_COLS <- c("#f5a623","#2780e3","#3fb618")
  
  observeEvent(input$compare_portfolios, {
    ms  <- get_mu_sigma()
    sel <- c(input$cmp_port1 %||% "none",
             input$cmp_port2 %||% "none",
             input$cmp_port3 %||% "none")
    
    # Resolve weights for each active slot
    ports <- list()
    for (i in seq_along(sel)) {
      if (sel[i] == "none") next
      w <- tryCatch(resolve_cmp_slot(sel[i], i), error=function(e) NULL)
      if (!is.null(w)) ports[[cmp_label(sel[i])]] <- w
    }
    
    if (length(ports) == 0) {
      showNotification("No valid portfolios selected for comparison.", type="warning"); return()
    }
    
    # ── Stats cards ──────────────────────────────────────────────
    output$compare_stats_ui <- renderUI({
      col_w <- floor(12 / length(ports))
      fluidRow(lapply(seq_along(ports), function(i) {
        nm <- names(ports)[i]; w <- ports[[i]]
        s  <- port_stats(w, ms$mu, ms$Sigma, input$rf_rate/100)
        set.seed(42)
        sims <- as.numeric(MASS::mvrnorm(5000,mu=ms$mu,Sigma=ms$Sigma) %*% w)
        q    <- quantile(sims, 0.05); cvar <- -mean(sims[sims<=q])*100
        acc  <- CMP_COLS[i]
        column(col_w, div(
          style=paste0("background:#1a1f2e;border:1px solid ",acc,
                       ";border-radius:4px;padding:14px 18px;margin-bottom:6px;"),
          tags$div(style=paste0("color:",acc,";font-family:'Source Sans 3';font-size:13px;",
                                "font-weight:600;margin-bottom:10px;border-bottom:1px solid ",
                                acc,";padding-bottom:6px;"), nm),
          tags$table(style="width:100%;",
                     tags$tr(
                       tags$td(style="color:#7a87a8;font-size:12px;padding:3px 0;","Return"),
                       tags$td(style="color:#d0d6e8;font-family:'IBM Plex Mono';font-size:13px;text-align:right;",
                               paste0(round(s$Return*100,2),"%"))),
                     tags$tr(
                       tags$td(style="color:#7a87a8;font-size:12px;padding:3px 0;","Volatility"),
                       tags$td(style="color:#d0d6e8;font-family:'IBM Plex Mono';font-size:13px;text-align:right;",
                               paste0(round(s$Vol,2),"%"))),
                     tags$tr(
                       tags$td(style="color:#7a87a8;font-size:12px;padding:3px 0;","Sharpe"),
                       tags$td(style=paste0("color:",acc,";font-family:'IBM Plex Mono';font-size:13px;text-align:right;"),
                               round(s$Sharpe %||% NA, 3))),
                     tags$tr(
                       tags$td(style="color:#7a87a8;font-size:12px;padding:3px 0;","CVaR (5%)"),
                       tags$td(style="color:#ff0039;font-family:'IBM Plex Mono';font-size:13px;text-align:right;",
                               paste0(round(cvar,2),"%")))
          )
        ))
      }))
    })
    
    # ── Grouped weights bar chart ────────────────────────────────
    output$compare_weights_plot <- renderPlotly({
      assets <- rv$assets$Asset
      fig <- plot_ly()
      for (i in seq_along(ports)) {
        nm <- names(ports)[i]; w <- ports[[i]]
        fig <- fig %>% add_trace(
          x=assets, y=round(w*100,2), type="bar", name=nm,
          marker=list(color=CMP_COLS[i], opacity=0.85, line=list(color="#0f1117",width=0.5)),
          hovertemplate=paste0("<b>",nm,"</b><br>%{x}: %{y:.2f}%<extra></extra>"))
      }
      fig %>% layout(barmode="group") %>%
        plot_dark("","Weight (%)", legend=TRUE)
    })
    
    # ── Side-by-side weights & risk table ───────────────────────
    output$compare_table <- renderDT({
      assets <- rv$assets$Asset
      df <- data.frame(Asset=assets, Class=rv$assets$Class, stringsAsFactors=FALSE)
      for (nm in names(ports)) {
        w  <- ports[[nm]]
        rc <- risk_contributions(w, ms$Sigma)
        df[[paste0(nm," Wt (%)")]]  <- round(w*100, 2)
        df[[paste0(nm," RC (%)")]]  <- round(rc$RC_pct, 2)
      }
      datatable(df, rownames=FALSE, class="compact stripe",
                options=list(dom="t", paging=FALSE, scrollX=TRUE))
    })
    
    showNotification(paste0("Comparing ", length(ports), " portfolios."), type="message")
  })
  
  ## 6. Black-Litterman ------------------------------------------------------
  
  output$bl_mkt_weights_ui <- renderUI({
    n <- nrow(rv$assets); eq <- round(100/n,1)
    lapply(seq_len(n), function(i) fluidRow(
      column(7,tags$label(style="color:#7a87a8;font-size:12px;line-height:34px;",rv$assets$Asset[i])),
      column(5,numericInput(paste0("bl_mkt_",i),NULL,value=eq,min=0,max=100,step=0.5))))
  })
  
  output$bl_views_ui <- renderUI({
    if (rv$bl_view_count==0) return(p(class="hint-text","No views added yet. Click Add View to begin."))
    n <- nrow(rv$assets)
    lapply(seq_len(rv$bl_view_count), function(i) {
      wellPanel(class="well-sm",
                tags$b(style="color:#d0d6e8;font-size:13px;font-weight:600;",paste("View",i)),
                p(class="hint-text","Asset weights: +long, -short, 0=excluded."),
                fluidRow(lapply(seq_len(n), function(j)
                  column(max(1,floor(12/n)),numericInput(paste0("bl_v_",i,"_w_",j),rv$assets$Asset[j],value=0,step=0.5)))),
                fluidRow(
                  column(6,numericInput(paste0("bl_v_",i,"_q"),"View Return (% p.a.)",value=1,step=0.5)),
                  column(6,numericInput(paste0("bl_v_",i,"_conf"),"Confidence (0-1)",value=0.5,min=0.01,max=1,step=0.05))))
    })
  })
  observeEvent(input$bl_add_view,    { rv$bl_view_count <- rv$bl_view_count+1 })
  observeEvent(input$bl_remove_view, { if (rv$bl_view_count>0) rv$bl_view_count <- rv$bl_view_count-1 })
  
  observeEvent(input$bl_run, {
    req(rv$bl_view_count>0); ms <- get_mu_sigma(); n <- nrow(rv$assets)
    w_mkt <- sapply(seq_len(n),function(i) input[[paste0("bl_mkt_",i)]]%||%(100/n))
    w_mkt <- w_mkt/sum(w_mkt)
    mu_eq <- bl_implied_returns(ms$Sigma,w_mkt,delta=input$bl_delta,rf=input$rf_rate/100)
    k <- rv$bl_view_count
    P <- matrix(0,k,n); Q <- numeric(k); od <- numeric(k)
    for (i in seq_len(k)) {
      for (j in seq_len(n)) P[i,j] <- input[[paste0("bl_v_",i,"_w_",j)]]%||%0
      Q[i]  <- (input[[paste0("bl_v_",i,"_q")]]%||%0)/100
      conf  <- (input[[paste0("bl_v_",i,"_conf")]]%||%0.5)
      od[i] <- ((1-conf)/conf)*input$bl_tau*as.numeric(P[i,]%*%ms$Sigma%*%P[i,])
    }
    od <- pmax(od,1e-10)
    bl <- tryCatch(bl_posterior(mu_eq,ms$Sigma,P,Q,omega=diag(od),tau=input$bl_tau),
                   error=function(e){showNotification(paste("BL error:",e$message),type="error");NULL})
    if (is.null(bl)) return()
    rv$bl_mu <- bl$mu; rv$bl_Sigma <- bl$Sigma
    wb <- get_weight_bounds(); gc <- get_all_constraints()
    fl <- generate_frontier(bl$mu,bl$Sigma,n_points=60,w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc)
    if (!is.null(fl)) rv$bl_weights <- max_sharpe_portfolio(fl,rf=input$rf_rate/100)$Weights
    ret_df <- data.frame(Asset=rv$assets$Asset,Equilibrium=round(mu_eq*100,2),BL=round(bl$mu*100,2))
    output$bl_returns_plot <- renderPlotly({
      plot_ly(ret_df) %>%
        add_trace(x=~Asset,y=~Equilibrium,type="bar",name="Equilibrium",
                  marker=list(color="#2780e3",opacity=0.65,line=list(color="#0f1117",width=1))) %>%
        add_trace(x=~Asset,y=~BL,type="bar",name="BL Posterior",
                  marker=list(color="#3fb618",line=list(color="#0f1117",width=1))) %>%
        layout(barmode="group") %>% plot_dark("","Expected Return (%)",legend=TRUE)
    })
    if (!is.null(rv$bl_weights)) {
      wd <- data.frame(Asset=rv$assets$Asset,Weight=round(rv$bl_weights*100,2))
      wd <- wd[wd$Weight>0.01,]
      output$bl_weights_plot <- renderPlotly({
        p <- plot_ly(wd,x=~Asset,y=~Weight,type="bar",
                     marker=list(color=PAL[seq_len(nrow(wd))],line=list(color="#0f1117",width=1)),
                     hovertemplate="%{x}: %{y:.2f}%<extra></extra>")
        plot_dark(p,ylab="Weight (%)")
      })
    }
    output$bl_table <- renderDT({
      out <- data.frame(Asset=rv$assets$Asset,Eq=round(mu_eq*100,2),BL=round(bl$mu*100,2),
                        Chg=round((bl$mu-mu_eq)*100,2))
      out$Dir <- ifelse(out$Chg>0,"▲",ifelse(out$Chg<0,"▼","—"))
      colnames(out) <- c("Asset","Equilib. Return (%)","BL Return (%)","Change (pp)","Dir.")
      datatable(out,rownames=FALSE,class="compact stripe",options=list(dom="t",paging=FALSE)) %>%
        formatStyle("Change (pp)",color=styleInterval(c(-0.001,0.001),c("#ff0039","#7a87a8","#3fb618")))
    })
    showNotification("Black-Litterman posterior computed.",type="message")
  })
  
  ## 7. Monte Carlo ----------------------------------------------------------
  
  observeEvent(input$run_mc, {
    ms <- get_mu_sigma(); gc <- get_all_constraints(); wb <- get_weight_bounds()
    w <- switch(input$mc_portfolio,
                "min_var"      = { req(rv$gmv_port);   rv$gmv_port$Weights },
                "max_sharpe"   = { req(rv$msr_port);   rv$msr_port$Weights },
                "max_div"      = max_div_portfolio(ms$Sigma,ms$vols,wb$w_min_vec,wb$w_max_vec,gc),
                "risk_parity"  = risk_parity_portfolio(ms$Sigma),
                "equal_weight" = equal_weight_portfolio(nrow(rv$assets)),
                "analysis"     = { req(rv$analysis_w); rv$analysis_w },
                "current"      = { w <- get_current_portfolio(); req(!is.null(w)); w },
                "bl_msr"       = { req(rv$bl_weights); rv$bl_weights })
    if (is.null(w)) { showNotification("Run Optimizer/Analysis first.",type="warning"); return() }
    sim_mu    <- if (input$mc_portfolio=="bl_msr"&&!is.null(rv$bl_mu))    rv$bl_mu    else ms$mu
    sim_Sigma <- if (input$mc_portfolio=="bl_msr"&&!is.null(rv$bl_Sigma)) rv$bl_Sigma else ms$Sigma
    hor <- input$mc_horizon; n_sim <- input$mc_sims
    init <- input$mc_initial; contrib <- input$mc_contrib; alpha <- input$mc_alpha/100
    withProgress(message="Simulating...",value=0.2,{
      res <- run_mc_sim(w,sim_mu,sim_Sigma,hor,n_sim,init,contrib); setProgress(1)
    })
    rv$mc_result <- res
    wl <- res$wealth; pr <- res$port_ret; terminal <- wl[,hor+1]
    q_a <- quantile(terminal,alpha); cvar_loss <- mean(init-terminal[terminal<=q_a])
    output$mc_vb_median <- renderValueBox({ valueBox(fmt_dollar(median(terminal)),"Median Terminal Value",color="yellow",icon=icon("coins")) })
    output$mc_vb_p10    <- renderValueBox({ valueBox(fmt_dollar(quantile(terminal,0.10)),"10th Pctl Value",color="blue",icon=icon("arrow-down")) })
    output$mc_vb_p90    <- renderValueBox({ valueBox(fmt_dollar(quantile(terminal,0.90)),"90th Pctl Value",color="green",icon=icon("arrow-up")) })
    output$mc_vb_cvar   <- renderValueBox({ valueBox(fmt_dollar(max(0,cvar_loss)),paste0("Exp. Loss | Bottom ",input$mc_alpha,"%"),color="red",icon=icon("exclamation-triangle")) })
    years <- 0:hor
    p5 <- apply(wl,2,quantile,.05); p25 <- apply(wl,2,quantile,.25)
    p50 <- apply(wl,2,quantile,.50); p75 <- apply(wl,2,quantile,.75)
    p95 <- apply(wl,2,quantile,.95)
    si  <- sample(n_sim,min(300,n_sim))
    output$mc_paths_plot <- renderPlotly({
      fig <- plot_ly()
      for (idx in si) fig <- fig %>% add_trace(x=years,y=wl[idx,],type="scatter",mode="lines",
                                               line=list(color="rgba(39,128,227,0.04)",width=0.5),showlegend=FALSE,hoverinfo="skip")
      fig %>%
        add_trace(x=years,y=p95,type="scatter",mode="lines",line=list(color="#3fb618",dash="dot",width=1.5),name="95th Pctl") %>%
        add_trace(x=years,y=p75,type="scatter",mode="lines",line=list(color="#3d9be9",dash="dash",width=1.5),name="75th Pctl") %>%
        add_trace(x=years,y=p50,type="scatter",mode="lines",line=list(color="#2780e3",width=2.5),name="Median") %>%
        add_trace(x=years,y=p25,type="scatter",mode="lines",line=list(color="#ff7518",dash="dash",width=1.5),name="25th Pctl") %>%
        add_trace(x=years,y=p5, type="scatter",mode="lines",line=list(color="#ff0039",dash="dot",width=1.5),name="5th Pctl") %>%
        plot_dark("Year","Portfolio Value ($)",legend=TRUE) %>%
        layout(yaxis=list(tickformat="$,.0f"))
    })
    output$mc_hist_plot <- renderPlotly({
      plot_ly(x=terminal/1e6,type="histogram",nbinsx=60,
              marker=list(color="#2780e3",opacity=0.75,line=list(color="#0f1117",width=0.5))) %>%
        add_segments(x=median(terminal)/1e6,xend=median(terminal)/1e6,y=0,yend=n_sim*.12,
                     line=list(color="#3fb618",width=2,dash="dash"),name="Median") %>%
        add_segments(x=q_a/1e6,xend=q_a/1e6,y=0,yend=n_sim*.12,
                     line=list(color="#ff0039",width=2,dash="dash"),name=paste0(input$mc_alpha,"th Pctl")) %>%
        plot_dark("Terminal Value ($M)","Frequency",legend=TRUE)
    })
    prl <- data.frame(Year=as.factor(rep(seq_len(hor),each=n_sim)),Return=as.vector(pr)*100)
    output$mc_box_plot <- renderPlotly({
      plot_ly(prl,x=~Year,y=~Return,type="box",fillcolor="rgba(39,128,227,0.2)",
              line=list(color="#2780e3"),marker=list(color="#2780e3",opacity=0.2,size=2)) %>%
        plot_dark("Year","Annual Return (%)")
    })
    # ── CAGR distribution ───────────────────────────────────────
    cagr <- apply(pr, 1, function(r) (prod(1 + r))^(1/hor) - 1) * 100
    cagr_p05 <- quantile(cagr, 0.05)
    cagr_p25 <- quantile(cagr, 0.25)
    cagr_p50 <- median(cagr)
    cagr_p75 <- quantile(cagr, 0.75)
    cagr_p95 <- quantile(cagr, 0.95)
    
    output$mc_cagr_plot <- renderPlotly({
      dens  <- density(cagr, n=512)
      y_max <- max(dens$y)
      plot_ly() %>%
        add_trace(x=dens$x[dens$x<=cagr_p05], y=dens$y[dens$x<=cagr_p05],
                  type="scatter", mode="lines", fill="tozeroy",
                  fillcolor="rgba(255,0,57,0.18)", line=list(color="transparent"),
                  name="Bottom 5%", showlegend=TRUE) %>%
        add_trace(x=dens$x[dens$x>=cagr_p95], y=dens$y[dens$x>=cagr_p95],
                  type="scatter", mode="lines", fill="tozeroy",
                  fillcolor="rgba(63,182,24,0.18)", line=list(color="transparent"),
                  name="Top 5%", showlegend=TRUE) %>%
        add_trace(x=dens$x[dens$x>=cagr_p25 & dens$x<=cagr_p75],
                  y=dens$y[dens$x>=cagr_p25 & dens$x<=cagr_p75],
                  type="scatter", mode="lines", fill="tozeroy",
                  fillcolor="rgba(39,128,227,0.18)", line=list(color="transparent"),
                  name="IQR (25-75%)", showlegend=TRUE) %>%
        add_trace(x=dens$x, y=dens$y, type="scatter", mode="lines",
                  line=list(color="#2780e3", width=2), name="Density", showlegend=FALSE) %>%
        add_segments(x=cagr_p05, xend=cagr_p05, y=0, yend=y_max*0.85,
                     line=list(color="#ff0039",width=1.5,dash="dot"), name="5th Pctl") %>%
        add_segments(x=cagr_p50, xend=cagr_p50, y=0, yend=y_max*0.95,
                     line=list(color="#3fb618",width=2,dash="dash"), name="Median") %>%
        add_segments(x=cagr_p95, xend=cagr_p95, y=0, yend=y_max*0.85,
                     line=list(color="#3fb618",width=1.5,dash="dot"), name="95th Pctl") %>%
        layout(annotations=list(
          list(x=cagr_p05, y=y_max*0.88, text=paste0(round(cagr_p05,1),"%"),
               showarrow=FALSE, font=list(color="#ff0039",family="IBM Plex Mono",size=11)),
          list(x=cagr_p50, y=y_max*0.98, text=paste0(round(cagr_p50,1),"%"),
               showarrow=FALSE, font=list(color="#3fb618",family="IBM Plex Mono",size=11)),
          list(x=cagr_p95, y=y_max*0.88, text=paste0(round(cagr_p95,1),"%"),
               showarrow=FALSE, font=list(color="#3fb618",family="IBM Plex Mono",size=11))
        )) %>%
        plot_dark("Annualised Return (% p.a.)","Density", legend=TRUE)
    })
    
    # ── Max drawdown and tail risk table ────────────────────────
    mdd <- apply(wl, 1, function(path) {
      peak <- cummax(path)
      max(ifelse(peak > 0, (peak - path) / peak, 0))
    }) * 100
    
    ci_levels <- c(0.90, 0.95, 0.99)
    risk_rows <- do.call(rbind, lapply(ci_levels, function(ci) {
      tail_p <- (1 - ci) / 2
      data.frame(
        CI           = paste0(ci*100, "%"),
        Best_Return  = paste0(round(quantile(cagr, 1 - tail_p), 2), "%"),
        Worst_Return = paste0(round(quantile(cagr, tail_p), 2), "%"),
        Worst_MDD    = paste0(round(quantile(mdd, ci), 2), "%"),
        stringsAsFactors=FALSE)
    }))
    colnames(risk_rows) <- c("Confidence Level","Best Annualised Return",
                             "Worst Annualised Return","Max Drawdown (Worst)")
    
    output$mc_risk_table <- renderDT({
      datatable(risk_rows, rownames=FALSE, class="compact stripe",
                options=list(dom="t", paging=FALSE, scrollX=TRUE)) %>%
        formatStyle("Best Annualised Return",  color="#3fb618") %>%
        formatStyle("Worst Annualised Return", color="#ff0039") %>%
        formatStyle("Max Drawdown (Worst)",    color="#c9a84c")
    })
  })
  
  ## 8. Comparison & Export --------------------------------------------------
  
  observeEvent(input$run_comparison, {
    req(rv$frontier); ms <- get_mu_sigma(); gc <- get_all_constraints(); wb <- get_weight_bounds()
    n <- nrow(rv$assets)
    w_list <- list(
      "Min Variance"        = rv$gmv_port$Weights,
      "Max Sharpe"          = rv$msr_port$Weights,
      "Max Diversification" = max_div_portfolio(ms$Sigma,ms$vols,wb$w_min_vec,wb$w_max_vec,gc),
      "Risk Parity"         = risk_parity_portfolio(ms$Sigma),
      "Equal Weight"        = equal_weight_portfolio(n),
      "Max Utility (l=3)"  = max_utility_portfolio(ms$mu,ms$Sigma,lambda=3,wb$w_min_vec,wb$w_max_vec,gc),
      "Min CVaR"            = withProgress(message="Min-CVaR...",value=0.5,{
        min_cvar_portfolio(ms$mu,ms$Sigma,w_min_vec=wb$w_min_vec,w_max_vec=wb$w_max_vec,group_constraints=gc)}),
      "Max Return"          = max_return_portfolio(ms$mu,wb$w_min_vec,wb$w_max_vec))
    w_cur <- get_current_portfolio()
    if (!is.null(w_cur))          w_list[["Current Portfolio"]] <- w_cur
    if (!is.null(rv$bl_weights))  w_list[["Black-Litterman"]]  <- rv$bl_weights
    if (!is.null(rv$analysis_w))  w_list[["Current Analysis"]] <- rv$analysis_w
    w_list <- w_list[!sapply(w_list,is.null)]
    set.seed(42); sc <- MASS::mvrnorm(4000,mu=ms$mu,Sigma=ms$Sigma)
    stats_df <- do.call(rbind,lapply(names(w_list),function(nm){
      w  <- w_list[[nm]]; s <- port_stats(w,ms$mu,ms$Sigma,input$rf_rate/100)
      pr <- as.numeric(sc %*% w); q <- quantile(pr,.05); cv <- -mean(pr[pr<=q])*100
      div <- as.numeric(t(w)%*%(ms$vols/100))/(s$Vol/100)
      data.frame(Portfolio=nm,Return=round(s$Return*100,2),Vol=round(s$Vol,2),
                 Sharpe=round(s$Sharpe%||%NA,3),CVaR_5=round(cv,2),DivRatio=round(div,3),stringsAsFactors=FALSE)
    }))
    wmat <- do.call(rbind,lapply(w_list,function(w) round(w*100,1)))
    rownames(wmat) <- names(w_list); colnames(wmat) <- rv$assets$Asset
    rv$comp_data <- list(stats=stats_df,weights=wmat,w_list=w_list)
  })
  
  output$comparison_stats_table <- renderDT({
    req(rv$comp_data); df <- rv$comp_data$stats
    colnames(df) <- c("Portfolio","Return (%)","Volatility (%)","Sharpe","CVaR 5% (%)","Div. Ratio")
    datatable(df,rownames=FALSE,class="compact stripe",options=list(dom="t",paging=FALSE)) %>%
      formatStyle("Return (%)",color="#2780e3") %>%
      formatStyle("Sharpe",color="#3fb618") %>%
      formatStyle("CVaR 5% (%)",color="#ff0039")
  })
  
  output$comparison_heatmap <- renderPlotly({
    req(rv$comp_data); wmat <- rv$comp_data$weights
    hm <- as.data.frame(wmat) %>% mutate(Portfolio=rownames(wmat)) %>%
      pivot_longer(-Portfolio,names_to="Asset",values_to="Weight")
    plot_ly(hm,x=~Asset,y=~Portfolio,z=~Weight,type="heatmap",
            colorscale=list(c(0,"#0f1117"),c(0.4,"#1a3a6b"),c(1,"#2780e3")),
            hovertemplate="%{y} -> %{x}<br>Weight: %{z:.1f}%<extra></extra>",
            text=~paste0(Weight,"%"),texttemplate="%{text}",showscale=TRUE) %>%
      layout(paper_bgcolor="#1a1f2e",plot_bgcolor="#0f1117",
             xaxis=list(tickfont=list(color="#7a87a8",family="Source Sans 3",size=12),side="bottom"),
             yaxis=list(tickfont=list(color="#d0d6e8",family="Source Sans 3",size=12)),
             margin=list(l=185,r=80,t=20,b=80))
  })
  
  ### Export Handlers ---------------------------------------------------------
  
  output$export_csv <- downloadHandler(
    filename=function() paste0("portfolio_comparison_",Sys.Date(),".csv"),
    content=function(file) {
      req(rv$comp_data); cd <- rv$comp_data
      wdf <- cbind(Portfolio=rownames(cd$weights),as.data.frame(cd$weights))
      write.csv(merge(wdf,cd$stats,by="Portfolio",all=TRUE),file,row.names=FALSE)
    })
  
  output$export_excel <- downloadHandler(
    filename=function() paste0("portfolio_optimizer_",Sys.Date(),".xlsx"),
    content=function(file) {
      req(rv$comp_data); cd <- rv$comp_data
      wdf    <- cbind(Portfolio=rownames(cd$weights),as.data.frame(cd$weights))
      cor_df <- cbind(Asset=rv$assets$Asset,as.data.frame(rv$cor_matrix))
      colnames(cor_df)[-1] <- rv$assets$Asset
      sheets <- list("Weights"=wdf,"Statistics"=cd$stats,"Assets"=rv$assets,"Correlations"=cor_df)
      if (!is.null(rv$bl_mu)&&!is.null(rv$bl_weights)) {
        ms <- get_mu_sigma(); wm <- rep(1/nrow(rv$assets),nrow(rv$assets))
        mq <- bl_implied_returns(ms$Sigma,wm,delta=input$bl_delta,rf=input$rf_rate/100)
        sheets[["Black-Litterman"]] <- data.frame(
          Asset=rv$assets$Asset,Equilibrium_Return=round(mq*100,4),
          BL_Posterior_Return=round(rv$bl_mu*100,4),BL_Weight=round(rv$bl_weights*100,2))
      }
      writexl::write_xlsx(sheets,path=file)
    })
  
  output$export_bl <- downloadHandler(
    filename=function() paste0("bl_analysis_",Sys.Date(),".csv"),
    content=function(file) {
      req(rv$bl_mu); ms <- get_mu_sigma()
      wm <- rep(1/nrow(rv$assets),nrow(rv$assets))
      mq <- bl_implied_returns(ms$Sigma,wm,delta=input$bl_delta,rf=input$rf_rate/100)
      df <- data.frame(Asset=rv$assets$Asset,
                       Equilibrium_Return=round(mq*100,4),
                       BL_Posterior_Return=round(rv$bl_mu*100,4),
                       Change_pp=round((rv$bl_mu-mq)*100,4),
                       BL_Weight=if (!is.null(rv$bl_weights)) round(rv$bl_weights*100,2) else NA)
      write.csv(df,file,row.names=FALSE)
    })
  
  output$export_mc <- downloadHandler(
    filename=function() paste0("monte_carlo_",Sys.Date(),".csv"),
    content=function(file) {
      req(rv$mc_result); wl <- rv$mc_result$wealth; hor <- ncol(wl)-1
      terminal <- wl[,hor+1]
      pd <- as.data.frame(wl); colnames(pd) <- paste0("Year_",0:hor)
      pd <- cbind(Simulation=seq_len(nrow(pd)),pd)
      write.csv(pd,file,row.names=FALSE)
    })
  
} # end server

shinyApp(ui, server)
