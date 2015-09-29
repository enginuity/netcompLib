\name{manOutputSettings}

\title{Output Settings}

\section{Purpose}{
The output adjustment for functions in this package (especially the simulation functions) have multiple parts to it, and has a common format for all functions. 
}

\section{Parameters}{
\itemize{
\item verbose [logical] :: Function prints progress if TRUE. In simple settings, this may be the only function parameter
\item vbset [vector[3]-int] :: The entries of this vector corresponds to: (Note that all negative values in this vector are treated as zeros)
\itemize{
\item vbset[1] - The number of functions and levels of within-function-calls that will produce output. eg. If this is set to 1, only this function (but none of the functions it calls) will produce output. 
\item vbset[2] - The number of functions and levels of within-function-calls that will produce detailed output. This should be at most the value of vbset[1]... 
\item vbset[3] - The number of filler characters. This is usually incremented by 2 on each function call (so within-function-calls have output that is indented)
}
}
}

\section{Programming Notes}{
When calling a function, vbset_new = vbset + c(-1, -1, 2). This says to reduce the number of levels of output in the called function, as well as increase the indentation. 
}
