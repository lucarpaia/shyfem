
c parameter for Venice Lagoon - venlag.h

c all dimensions give the maximum possible dimension
c actual used values may be lower

	integer nkndim			!total number of points
	integer neldim			!total number of elements
	integer nlvdim			!total number of vertical levels
	integer mbwdim			!maximum bandwidth
	integer ngrdim			!maximum grade of node
	integer nbcdim			!total number of open boundaries
	integer nrbdim			!total number of open boundary nodes
	integer nardim			!maximum area code
	integer nexdim			!total number of extra points
	integer nfxdim			!total number of flux points

	parameter (nkndim=6000)
	parameter (neldim=12000)
	parameter (nlvdim=1)
	parameter (mbwdim=90)
	parameter (ngrdim=11)
	parameter (nbcdim=15)
	parameter (nrbdim=200)
	parameter (nardim=12)
	parameter (nexdim=100)
	parameter (nfxdim=1000)

c do not change anything after this line

	integer nlidim			!dimension for side index
        parameter (nlidim=2*nkndim+neldim)

c old unused dimension declarations

	integer nfldim,nfddim
	integer ipcdim,ivcdim
	integer ndldim
	parameter (nfldim=1,nfddim=1)
	parameter (ipcdim=1,ivcdim=1)
	parameter (ndldim=1)

