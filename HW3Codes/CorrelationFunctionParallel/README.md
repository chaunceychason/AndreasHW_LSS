# READ ME FOR CORRELATION FUNCTION STUFF


So there is lots that I am not going to say because it is 3:30am and Im pretty tired. SO here is the nitty gritty I might forget:

The files labelled output***March23.txt have the results from the serial-C code (misnamed HW3ParallelCorrelationFunc1_1.c). 
	The serial code double counted in the loops but did NOT correct for this in the counts, thus they are off by *2.
The files labelled output***March23Para.txt have the results from the semiparallel-C code (misnamed HW3ParallelCorrelationFunc2_1.c). 
	This can be used for comparason. The code here took 4 minutes to fun in full. 

The files without a data are the most current fully parallel codes. These correct by the division and output the results. 



