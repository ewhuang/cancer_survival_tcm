library(ggplot2)

# Plot the SCC data points.
df <- data.frame(chi=c(4.921097, 6.700698, 10.79155, 10.22641, 5.812572,
    6.214136, 6.214136, 6.214136, 6.214136, 6.214136), num_dim=c(300, 400, 500,
    600, 700), meth=c('HEMnet', 'HEMnet', 'HEMnet', 'HEMnet', 'HEMnet', 
    'Baseline', 'Baseline', 'Baseline', 'Baseline', 'Baseline'))

pl <- ggplot(data=df, aes(x=num_dim, y=chi, group=meth)) + geom_line(aes(
    color=meth)) + labs(x='Number of dimensions', y='Chi-squared statistic',
colour='Method type') + ylim(0, 12) + xlim(300, 700) + theme(aspect.ratio=0.66,
    text=element_text(size=12), legend.justification=c(1,0),
    legend.position=c(1,0), legend.background = element_rect(fill=alpha(
        'white', 0.8)))

ggsave(filename='./results/scc_sensitivity.pdf', plot=pl, width=5, height=3.3)

# Plot the non-sq nsclc data points.
df <- data.frame(chi=c(5.824815, 7.800769, 8.449226, 6.95555, 1.170268,
    5.143644, 5.143644, 5.143644, 5.143644, 5.143644), num_dim=c(300, 400, 500,
    600, 700), meth=c('HEMnet', 'HEMnet', 'HEMnet', 'HEMnet', 'HEMnet',
    'Baseline', 'Baseline', 'Baseline', 'Baseline', 'Baseline'))

pl <- ggplot(data=df, aes(x=num_dim, y=chi, group=meth)) + geom_line(aes(
    color=meth)) + labs(x='Number of dimensions', y='Chi-squared statistic',
colour='Method type') + ylim(0, 12) + xlim(300, 700) + theme(aspect.ratio=0.66,
    text=element_text(size=12), legend.justification=c(1,1),
    legend.position=c(1,1), legend.background = element_rect(fill=alpha(
        'white', 0.8)))

ggsave(filename='./results/nonsq_sensitivity.pdf', plot=pl, width=5, height=3.3)