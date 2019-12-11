from matplotlib import pyplot

num_procs = [ 1 , 2 , 4 , 6 , 8 , 10 , 12 , 14 , 16 ]

times_arrays = []
speedups_arrays = []
labels = []

times_arrays.append( [ 15.06 , 7.66 , 3.94 , 2.67 , 2.06 , 1.85 , 2.19 , 2.26 , 2.38 ] ) , labels.append("")


for time_array in times_arrays:
    speedups_arrays.append( [ (time_array[0]/x) for x in time_array ] )


pyplot.xlabel("Processadores")
pyplot.ylabel("Speedup")
pyplot.xlim(0,20)
pyplot.ylim(0,20)
pyplot.xticks(range(0,20,1))
pyplot.yticks(range(0,20,1))
pyplot.plot( num_procs , num_procs , color=(0,0,0,0.3) )

for i in range(0,len(speedups_arrays)):
    speedup_array =  speedups_arrays[i]
    pyplot.plot( num_procs , speedup_array , alpha=0.6 , label=labels[i] )
    pyplot.scatter( num_procs , speedup_array , alpha=1.0 )

pyplot.legend()
pyplot.show()