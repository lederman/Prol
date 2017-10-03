# A Note on Eigendecomposition

## General comments

The preliminary version of this code uses Matlab's eigendecomposition for readability. This particular implementation of eigndecomposition is not ideal for the task for several obvious reasons
* It requires us to construct the full matrix, although the matrix is tridiagonal, and although matlab uses the tridiagonal structure internally.
* We do not need all the eigenvectors, only the smaller ones in magnitude. 
There is another, more subtle reason, discussed next.

## Not all equals are equal

Generally speaking, the eigenvectors in an eigendecomposition are computed up to some accuracy &#949;, 
in the sense that if v is the correct normalized eigenvector and \||v-u\||<&#949; (in the natural l<sup>2</sup> norm), then the vector u is considered to be an accurate solution. 
In a well-conditioned problem, &#949; is "machine precision" at best (roughly 10<sup>-16</sup>). 
Therefore, if some entry of the eigenvector v is much smaller than &#949;, we might as well set it to zero; 
Suppose that the first element v(0) in the the eigenvector v is 10<sup>-30</sup>. Obviously, we can set u(0)=0 and still have u being a numerically accurate eigenvector. 

However, it has been shown in \[1\] that in some cases is it possible to achieve much higher accuracy in the computation of certain elements of the eigenvectors. 
This fact turns out to be very useful in several different computations associated with the eigenvalues of prolates. 
While many eigendecomposition algorithms obtain accurate eignvectors, not all algorithms retrieve these small elements to the desired accuracy, thus adversely effecting our computation. 
In this sense, many algorithms produce numerically "equal" eigenvectors (up to machine precision and sign/phase), but not all have this special property. 
In the preliminary version of this code, we use matlab's eigendecomposition in order to make the code more readable. We use an inverse power method step obtain an eigenvector with the desired properties (in some regimes) from the original eigenvector computed by matlab. 
In future versions we will replace this operation with an independent eigendecomposition. 


\[1\] Osipov, Andrei. "Small coordinates of eigenvectors of certain symmetric tridiagonal matrices: numerical evaluation and error analysis." (2014).

