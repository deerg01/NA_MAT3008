#gnuplot ./Homework#4/Code/genHisto.plt 

# 전체 그림의 크기 설정
set terminal png size 1200,600 enhanced
set output 'histograms.png'

# 기본 스타일 설정
set style fill solid 0.5
set boxwidth 0.4 relative
set grid ytics

# Multipanels 설정
set multiplot layout 2, 4 title "Distributions Histograms"

# 균일 분포 히스토그램
do for [sample in '100 1000 10000 100000'] {
    set title sprintf('Uniform Distribution (%s Samples)', sample)
    set xlabel 'x axis'
    set ylabel 'Probability Density'
    
    # 히스토그램 그리기
    plot sprintf('u_%s.txt', sample) using 1:(1.0) smooth freq with boxes lc rgb 'blue' title 'Uniform Distribution'
}

# 가우시안 분포 히스토그램
do for [sample in '100 1000 10000 100000'] {
    set title sprintf('Gaussian Distribution (%s Samples)', sample)
    set xlabel 'x axis'
    set ylabel 'Probability Density'
    
    # 히스토그램 그리기
    plot sprintf('g_%s.txt', sample) using 1:(1.0) smooth freq with boxes lc rgb 'red' title 'Gaussian Distribution'
}

unset multiplot
