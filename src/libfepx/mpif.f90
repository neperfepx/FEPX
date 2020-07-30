!@@
! MPI/Pro and related software.
! Copyright (C) 1997-2004 MPI Software Technology, Inc.  All rights reserved.
! 
! DISCLAIMER OF LIABILITY:
! 
! THIS SOFTWARE IS PROVIDED BY MPI SOFTWARE TECHNOLOGY, INC.
! AND ITS LICENSORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES,
! INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
! AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL
! MPI SOFTWARE TECHNOLOGY, INC. NOR ITS LICENSORS BE LIABLE FOR ANY DIRECT,
! INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES
! (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
! SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
! HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT,
! STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING
! IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF MPI SOFTWARE
! TECHNOLOGY, INC. HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGES.
!
! Use of this software in source code form is subject to written source 
! code to a written source code licensing agreement ("SCLA") between 
! MPI Software Technology, Inc. and the LICENSEE and sub-LICENSEE, 
! which define the rights and responsibilities of both LICENSEE and 
! sub-LICENSEE, and limits MPI Software Technology's liability. 
! This copyright notice must not be removed on copies and derivative 
! works made under LICENSE or sub-LICENSE. No use of the source code of 
! this product without a source code license agreement is permitted. 
! All derivative works become the property of MPI Software Technology, Inc.
! 
! Other copyrights may apply to this software, and are noted where
! appropriate.
! 
! MPI Software Technology, Inc.
! 110 12th Street North, Suite D103
! Birmingham, AL 35203
! TEL: +1 (205)314-3471, FAX: +1 (205) 314-3474
! http://www.mpi-softtech.com
!@!

Module MPI

Integer, parameter :: MPI_SUCCESS		= 0
Integer, parameter :: MPI_ERR_BUFFER		= 1
Integer, parameter :: MPI_ERR_COUNT		= 2
Integer, parameter :: MPI_ERR_TYPE		= 3
Integer, parameter :: MPI_ERR_TAG		= 4
Integer, parameter :: MPI_ERR_COMM		= 5
Integer, parameter :: MPI_ERR_RANK		= 6
Integer, parameter :: MPI_ERR_ROOT		= 7
Integer, parameter :: MPI_ERR_GROUP		= 8
Integer, parameter :: MPI_ERR_OP		= 9
Integer, parameter :: MPI_ERR_TOPOLOGY		= 10
Integer, parameter :: MPI_ERR_DIMS		= 11
Integer, parameter :: MPI_ERR_ARG		= 12
Integer, parameter :: MPI_ERR_UNKNOWN		= 13
Integer, parameter :: MPI_ERR_TRUNCATE		= 14
Integer, parameter :: MPI_ERR_OTHER		= 15
Integer, parameter :: MPI_ERR_INTERN		= 16
Integer, parameter :: MPI_ERR_IN_STATUS 	= 17
Integer, parameter :: MPI_ERR_PENDING		= 18
Integer, parameter :: MPI_ERR_REQUEST		= 19
Integer, parameter :: MPI_ERR_FILE		= 20
Integer, parameter :: MPI_ERR_UNSUPPORTED_OPERATION = 21
Integer, parameter :: MPI_ERR_IO		= 22
Integer, parameter :: MPI_ERR_AMODE		= 23
Integer, parameter :: MPI_ERR_UNSUPPORTED_DATAREP = 24
Integer, parameter :: MPI_ERR_NO_MEM		= 25
Integer, parameter :: MPI_ERR_BASE		= 26
Integer, parameter :: MPI_ERR_NOT_INITIALIZED 	= 27
Integer, parameter :: MPI_ERR_FINALIZED	  	= 28
Integer, parameter :: MPI_ERR_COMM_MEDIUM 	= 29
Integer, parameter :: MPI_ERR_NODE_DOWN 	= 30
Integer, parameter :: MPI_ERR_CANCELLED		= 31
Integer, parameter :: MPI_ERR_LASTCODE		= 32

Integer, parameter :: MPI_CART = 1
Integer, parameter :: MPI_GRAPH = 2
Integer, parameter :: MPI_BSEND_OVERHEAD = 32
Integer, parameter :: MPI_SOURCE = 2
Integer, parameter :: MPI_TAG = 3
Integer, parameter :: MPI_ERROR = 4
Integer, parameter :: MPI_STATUS_SIZE = 4
Integer, parameter :: MPI_MAX_PROCESSOR_NAME = 80
Integer, parameter :: MPI_MAX_ERROR_STRING = 80
Integer, parameter :: MPI_MAX_NAME_STRING = 80
Integer, parameter :: MPI_MAX_OBJECT_NAME = 128

Integer, parameter :: MPI_ERRHANDLER_NULL = 0
Integer, parameter :: MPI_ERRORS_ARE_FATAL = 1
Integer, parameter :: MPI_ERRORS_RETURN = 2

Integer, parameter :: MPI_KEYVAL_INVALID = 0
Integer, parameter :: MPI_TAG_UB = 1
Integer, parameter :: MPI_HOST = 2
Integer, parameter :: MPI_IO = 3
Integer, parameter :: MPI_WTIME_IS_GLOBAL = 4
Integer, parameter :: MPI_NULL_COPY_FN = 1
Integer, parameter :: MPI_DUP_FN = 2
Integer, parameter :: MPI_NULL_DELETE_FN = 3

Integer, parameter :: MPI_REQUEST_NULL = 0
Integer, parameter :: MPI_INFO_NULL = 0
Integer, parameter :: MPI_COMM_NULL = 0
Integer, parameter :: MPI_COMM_WORLD = 1
Integer, parameter :: MPI_COMM_SELF = 2
Integer, parameter :: MPI_GROUP_NULL = 0
Integer, parameter :: MPI_GROUP_EMPTY = 1

Integer, parameter :: MPI_UNEQUAL=0
Integer, parameter :: MPI_IDENT=1
Integer, parameter :: MPI_CONGRUENT=2
Integer, parameter :: MPI_SIMILAR=3

Integer, parameter :: MPI_DATATYPE_NULL=0
Integer, parameter :: MPI_PACKED = 14
Integer, parameter :: MPI_LB = 15
Integer, parameter :: MPI_UB = 16
Integer, parameter :: MPI_CHARACTER = 23
Integer, parameter :: MPI_LOGICAL = 24
Integer, parameter :: MPI_INTEGER = 25
Integer, parameter :: MPI_REAL = 26
Integer, parameter :: MPI_DOUBLE_PRECISION = 27
Integer, parameter :: MPI_COMPLEX = 28
Integer, parameter :: MPI_DOUBLE_COMPLEX = 29
Integer, parameter :: MPI_2INTEGER = 30
Integer, parameter :: MPI_2REAL = 31
Integer, parameter :: MPI_2DOUBLE_PRECISION = 32
Integer, parameter :: MPI_2COMPLEX = 33
Integer, parameter :: MPI_2DOUBLE_COMPLEX = 34
Integer, parameter :: MPI_INTEGER1 = 1
Integer, parameter :: MPI_INTEGER2 = 4
Integer, parameter :: MPI_INTEGER4 = 25
Integer, parameter :: MPI_BYTE = 3
Integer, parameter :: MPI_REAL4 = 26
Integer, parameter :: MPI_REAL8 = 27

Integer, parameter :: MPI_OP_NULL=0
Integer, parameter :: MPI_MAX=1
Integer, parameter :: MPI_MIN=2
Integer, parameter :: MPI_SUM=3
Integer, parameter :: MPI_PROD=4
Integer, parameter :: MPI_LAND=5
Integer, parameter :: MPI_BAND=6
Integer, parameter :: MPI_LOR=7
Integer, parameter :: MPI_BOR=8
Integer, parameter :: MPI_LXOR=9
Integer, parameter :: MPI_BXOR=10
Integer, parameter :: MPI_MAXLOC=11
Integer, parameter :: MPI_MINLOC=12

Integer, parameter :: MPI_PLEVEL_DISABLED = 0
Integer, parameter :: MPI_PLEVEL_NORMAL = 1
Integer, parameter :: MPI_PLEVEL_BUFFERS_FLUSHED = 2

Integer, parameter :: MPI_ANY_TAG = (-1)
Integer, parameter :: MPI_ANY_SOURCE = (-2)
Integer, parameter :: MPI_PROC_NULL	= (-1)
Integer, parameter :: MPI_UNDEFINED	= (-32766)
Integer, parameter :: MPI_VERSION = 1
Integer, parameter :: MPI_SUBVERSION = 2

Common /MPIPRO/ MPI_BOTTOM
Integer MPI_BOTTOM
Save /MPIPRO/

double precision :: MPI_WTIME, MPI_WTICK
external MPI_WTIME, MPI_WTICK

End Module MPI

