## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo=TRUE)
options(knitr.table.format = "html")

## ----how_to_install, eval=FALSE-----------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("williazo/mvGPS")

## ----dag_draw, echo=FALSE, fig.align="center", out.width="50%", cache=TRUE----
require(dagitty)
require(ggdag)
dag_coords <- data.frame(x=c(5.5, 5.5, 10.5, 8, 8, 10.5), 
                         y=c(4, 7, 7, 6, 3, 4), 
                         name=c("C1", "D1", "D2", "C2", "Y", "C3"))
dag_saturated <- dagify(Y~D1+D2+C1+C2+C3,
                        D1~C1+C2, D2~C2+C3,
                        D1~~D2, outcome="Y", coords=dag_coords)
sim_dag <- dag_saturated %>% 
    tidy_dagitty(layout = "nicely", seed = 12345) %>%
    dplyr::arrange(name) %>% 
    ggplot(aes(x = x, y = y, xend = xend, yend = yend)) +
    geom_dag_edges() +
    geom_dag_text(parse = TRUE, 
                  label = c(bquote(C[1]), bquote(C[2]), bquote(C[3]), bquote(D[1]), bquote(D[2]), "Y"), size=10, 
                  color=c("black")) +
    theme_dag()
sim_dag

## ----gen_exposure, message=FALSE, warning=FALSE-------------------------------
require(mvGPS)
sim_dt <- gen_D(method="u", n=200, rho_cond=0.2, s_d1_cond=2, s_d2_cond=2, k=3, 
                C_mu=rep(0, 3), C_cov=0.1, C_var=1, 
                d1_beta=c(0.5, 1, 0), d2_beta=c(0, 0.3, 0.75), seed=06112020)
D <- sim_dt$D
C <- sim_dt$C

## ----marg_param, echo=FALSE---------------------------------------------------
rho <- sim_dt$rho

## ----y_gen--------------------------------------------------------------------
alpha <- c(0.75, 1, 0.6, 1, 1)
sd_Y <- 2
X <- cbind(C, D)
Y <- X%*%alpha + rnorm(200, sd=sd_Y)

## ----mvGPS_w, message=FALSE---------------------------------------------------
require(mvGPS)
out_mvGPS <- mvGPS(D=D, C=list(C[, 1:2], C[, 2:3]))
w <- out_mvGPS$w

## ----balance_summary, warning=FALSE, message=FALSE----------------------------
require(knitr)
bal_results <- bal(model_list=c("mvGPS", "entropy", "CBPS", "PS", "GBM"), D, C=list(C[, 1:2], C[, 2:3]))
bal_summary <- bal_results$bal_metrics 
#contains overall summary statistics with respect to balance
bal_summary <-data.frame(bal_summary, ESS=c(bal_results$ess, nrow(D)))
#adding in ESS with last value representing the unweighted case
bal_summary <- bal_summary[order(bal_summary$max_cor), ]

kable(bal_summary[, c("euc_dist", "max_cor", "avg_cor", "ESS", "method")], 
      digits=4, row.names=FALSE, 
      col.names=c("Euc. Distance", "Max. Abs. Corr.", 
                  "Avg. Abs. Corr.", "ESS", "Method"))

## ----bias_est, warning=FALSE, message=FALSE-----------------------------------
dt <- data.frame(Y, D)
mvGPS_mod <- lm(Y ~ D1 + D2, weights=w, data=dt)
mvGPS_hat <- coef(mvGPS_mod)[c("D1", "D2")]

unadj_hat <- coef(lm(Y ~ D1 + D2, data=dt))[c("D1", "D2")]

bias_tbl <- cbind(truth=c(1, 1), unadj=unadj_hat, mvGPS_hat)
kable(bias_tbl, digits=2, row.names=TRUE, 
             col.names=c("Truth", "Unadjusted", "mvGPS"))

## ----chull_plot, warning=FALSE, message=FALSE, fig.align="center", out.width="50%", cache=TRUE----
require(sp)
chull_D <- hull_sample(D) 
#generate convex hull of exposure
chull_D_trim <- hull_sample(D, trim_hull=TRUE, trim_quantile=0.95)
#generate trimmed convex hull

bbox_grid <- sp::bbox(chull_D$hpts_vs) #bounding box over convex hull
bbox_df <- data.frame(D1=c(bbox_grid[1, 1], bbox_grid[1, 2], 
                           bbox_grid[1, 2], bbox_grid[1, 1]), 
                      D2=c(bbox_grid[2, 1], bbox_grid[2, 1], 
                           bbox_grid[2, 2], bbox_grid[2, 2]))
bbox_grid_trim <- sp::bbox(chull_D_trim$hpts_vs) #bounding box over trimmed convex hull
bbox_df_trim <- data.frame(D1=c(bbox_grid_trim[1, 1], bbox_grid_trim[1, 2],
                                bbox_grid_trim[1, 2], bbox_grid_trim[1, 1]), 
                           D2=c(bbox_grid_trim[2, 1], bbox_grid_trim[2, 1], 
                                bbox_grid_trim[2, 2], bbox_grid_trim[2, 2]))
chull_plot <- ggplot(data=data.frame(D), aes(x=D1, D2))+
    geom_point()+
    geom_polygon(data=data.frame(chull_D$hpts_vs), color="indianred4", fill=NA)+
    geom_polygon(data=data.frame(chull_D_trim$hpts_vs), color="indianred1", fill=NA, alpha=0.4)+
    geom_polygon(data=bbox_df, color="dodgerblue4", fill=NA)+
    geom_polygon(data=bbox_df_trim, color="dodgerblue1", fill=NA, alpha=0.4)+
    xlab("D1")+
    ylab("D2")+
    theme_bw()
chull_plot

