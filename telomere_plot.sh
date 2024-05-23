tidk search -s $2 -o $3 -d . $1
tidk plot -o $3_plot -t {$3}_telomeric_repeat_windows.tsv
convert {$3}_plot.svg {$3}_plot.jpg
