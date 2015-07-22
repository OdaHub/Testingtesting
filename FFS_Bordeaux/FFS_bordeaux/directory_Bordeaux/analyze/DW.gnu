set tit "DoubleWell(x)"                                                         
set xran [-2:2]                                                                 
set yran [-1:.10]                                                               
set xlab "x"                                                                    
set ylab " "                                                                    
 k4=   1.0000000000000000                                                       
 k2=  -2.0000000000000000                                                       
                                                                                
                                                                                
V(x)=k4*x**4+k2*x**2                                                            
plot V(x)                                                                       
                                                                                
                                                                                
 set ter pos  enhance col sol
 set ou "DW.eps"
 repl
 set ter x11
