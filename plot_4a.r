library(tidyverse)
library(ggplot2)
library(dplyr)
library(patchwork)
library(ggpubr)
library(gridExtra)


ac_da_all_mean <- read.csv('./files/genetic_single_performance_top100.csv')
method_order <- c('CPA', 'biolord', 'scouter', 'AttentionPert','baseMLP', "scELMo" ,"GeneCompass",'scGPT',
                  'scFoundation', 'GEARS','GenePert', 'linearModel', 'trainMean', 'controlMean', "baseReg" )
ac_da_all_mean$method <- factor(ac_da_all_mean$method, levels = rev(method_order))

ac_da_all2 <- ac_da_all_mean[, c('method', 'ac_score', 'Rank_all_final', 'DataSet')]
ac_da_all2$plot_name = 'all'

ac_da_mse <- ac_da_all_mean[, c('method', 'mse_score', 'Rank_mse_fianl', 'DataSet')]
ac_da_mse$plot_name = 'mse'
ac_da_pcc <- ac_da_all_mean[, c('method', 'cor', 'Rank_pcc_fianl', 'DataSet')]
ac_da_pcc$plot_name = 'PCC'
ac_da_edis <- ac_da_all_mean[, c('method', 'edistance_score', 'Rank_edistance_fianl', 'DataSet')]
ac_da_edis$plot_name = 'edis'

ac_da_deg <- ac_da_all_mean[, c('method', 'DEG_score', 'Rank_deg_fianl', 'DataSet')]
ac_da_deg$plot_name = 'DEG'
ac_da_was <- ac_da_all_mean[, c('method', 'was_score', 'Rank_was_final', 'DataSet')]
ac_da_was$plot_name = 'was'
ac_da_sym <- ac_da_all_mean[, c('method', 'sym_score', 'Rank_sym_final', 'DataSet')]
ac_da_sym$plot_name = 'sym'

p_list <- list()


df_list <- list(ac_da_all2, ac_da_mse, ac_da_pcc, ac_da_edis, ac_da_deg, ac_da_was, ac_da_sym)
for (df in df_list) {
  pname <- unique(df$plot_name)
  colnames(df) <- c('method', 'score', 'Rank', 'DataSet', 'plot_name')
  
  df <- df %>%
    group_by(method) %>%
    mutate(mean_rank = mean(Rank, na.rm = TRUE)) %>% 
    ungroup() 
  
  df <- df %>%
    group_by(DataSet) %>%
    mutate(final_rank = rank(mean_rank, ties.method = "min")) %>%
    data.frame()
  
  p = ggplot(data = df, aes(y = method, x = score)) +
    stat_summary(aes(fill = final_rank), 
                 fun = "mean", 
                 geom = "bar", 
                 width = 0.8, 
                 color = 'black', 
                 size = 0.2) +
    geom_point(position = position_jitter(width = 0, height = 0.2), 
               size = 0.15, alpha = 0.5) +
    stat_summary(fun.data = 'mean_se', 
                 geom = "errorbar", 
                 colour = "black", 
                 width = 0.3) +
    
    scale_fill_gradient(low = "white", high = "navyblue") + 
    geom_text(data = filter(df, final_rank <= 3),
              aes(x = 0.5, y = method,label = final_rank)) +
    labs(x = 'Score', y = 'Methods', title = pname) +
    theme(panel.grid.major=element_line(colour=NA),
          panel.background = element_rect(fill = "transparent",colour = 'black' ),
          panel.border = element_blank(), 
          axis.ticks = element_blank(), 
          plot.background = element_rect(fill = "white"),
          panel.grid.minor = element_blank(),
          plot.title = element_text(hjust = 0.5, size = 20, color = "black"),
          legend.text = element_text(size = 16, color = "black"),
          legend.title = element_text(size = 16, color = "black"),
          #axis.text = element_blank(),
          axis.title = element_blank())+
    xlim(0, 1) +
    theme(legend.position = "none")
  
  p_list[[pname]] <- p
}

pall = grid.arrange(grobs = p_list, ncol = 7)
