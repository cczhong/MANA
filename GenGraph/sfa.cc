#include "sfa.h"

using namespace std;

SFA::SFA()
{
    init();
}

SFA::SFA( char **s, int n )
{
    init();
    seqs = s;
    nreads = n;
}

SFA::~SFA()
{
    clear();
}

void SFA::init()
{
	concat  = NULL;
	SA      = NULL;
	LCP     = NULL;
	mLCP    = NULL;
	seqs    = NULL;
	nreads  = 0;
    verbose = false;
}

void SFA::clear()
{
    purgeSA();
    purgeDoc();
    purgeLCP();
    purgeMLCP();
}

void SFA::purgeSA()
{
    if ( SA != NULL ) delete[] SA;
    SA = NULL;
}

void SFA::purgeDoc()
{
    if ( concat != NULL ) delete[] concat;
    concat = NULL;
}

void SFA::purgeLCP()
{
    if ( LCP != NULL ) delete[] LCP;
    LCP = NULL;
}

void SFA::purgeMLCP()
{
    mLCP = NULL;
}


void  SFA::setStorageSize()
{
    //size_t max_seq = std::numeric_limits<LCPTYPE>::max();
    //size_t max_doc = std::numeric_limits<SFAIDXTYPE>::max();

    size = 0;
    for ( int i = 0; i < nreads; i++ ) {
        int len = strlen(seqs[i]);
        /*
        if ( len >= max_seq ) {
            printf("Long sequence given:%zu. Maximum sequence length:%zu\n", len, max_seq);
            exit(EXIT_FAILURE);
        }
        */
        size += len;
        size++; // terminater 
    }
    
    /*
    if ( size > max_doc ) {
        printf("Lengthy input:%zu. Maximum document size:%zu\n", size, max_doc);
        exit(EXIT_FAILURE);
    }
    */
}


void SFA::concatenateAllSeqs()
{
    SFAIDXTYPE cpos = 0;
    for ( int i = 0; i < nreads; i++ ) {
        int len = strlen(seqs[i]);
        memcpy( &concat[cpos], (SFACHARTYPE *)seqs[i], len );
        cpos += len;
        concat[cpos] = DELIMITER;
        cpos++;
    }
    
    assert(size == cpos);
    concat[size] = '\0';
}

void SFA::buildSFA()
{
    setStorageSize( );

    concat = new SFACHARTYPE[size+1];
    if ( concat == NULL ) {
        std::cerr << "Memory allocation error\n";
        exit(EXIT_FAILURE);
    }
	
    concatenateAllSeqs( );
    
    // DEBUG
    //cout << "Concatenated sequence: " << concat << endl;
    
    assert(strlen((char*) concat) > 0);
    
    SA = new SFAIDXTYPE[size];
    if ( concat == NULL ) {
        std::cerr << "Memory allocation error\n";
        exit(EXIT_FAILURE);
    }
    
    divsufsort( concat, SA, size );
// #if BIGINT == 0
//    divsufsort( concat, SA, size );
// #else
//     divsufsort64( concat, SA, size );
// #endif
}

// Simple Kasai
void SFA::buildLCP()
{
    assert(SA != NULL);
    SFAIDXTYPE *rank = new SFAIDXTYPE[size];
    for (SFAIDXTYPE i = 0; i < size; i++ )
        rank[SA[i]] = i;

    LCP = new LCPTYPE[size];
    LCP[0] = 0;
    SFAIDXTYPE h = 0;
    for (SFAIDXTYPE i = 0 ; i < size; i++ ) {
        if ( rank[i] > 0 ) {
            SFAIDXTYPE j = SA[rank[i]-1];
            while( concat[i+h] != DELIMITER && concat[i+h] == concat[j+h] ) {
                h++;
                assert(i+h<=size && j+h<=size);
            }
            LCP[rank[i]] = h;
            //LCP[rank[i]-1] = h;
            if ( h > 0 ) h--;
        }
    }
    //LCP[size-1] = 0;
    delete[] rank;
}


void SFA::buildMLCP()
{
    ia.build( LCP, size );

    SFAIDXTYPE n = ia.getSize();
    assert( n == size-1 );

    mLCP = ia.getIntervals();
}


void SFA::writeSFA( const char *filename )
{
    fio::write<SFAIDXTYPE>( SA, size, filename );
}

void SFA::writeDoc( const char *filename )
{
    fio::write<SFACHARTYPE>( concat, size, filename );
}

void SFA::writeLCP( const char *filename )
{
    if ( LCP == NULL ) return;
    fio::write<LCPTYPE>( LCP, size, filename );
}

void SFA::writeMLCP( const char *filename )
{
    if ( LCP == NULL ) return;
    ia.dump(filename);
}

void SFA::dump( const char *sfa_file,
                const char *doc_file,
                const char *lcp_file,
                const char *mcp_file )
{
    //writeSFA( sfa_file );
    //writeDoc( doc_file );
    writeLCP( lcp_file );
    writeMLCP( mcp_file );
}

void SFA::printSFA() 
{
    for (SFAIDXTYPE i = 0; i < size; i++ ) 
        std::cout << SA[i] << " ";
    std::cout << "\n";
}

void SFA::printDoc()
{
    std::cout << concat << "\n";
}

void SFA::printLCP()
{
    for (SFAIDXTYPE i = 0; i < size; i++ ) 
        std::cout << (unsigned)LCP[i] << " ";
    std::cout << "\n";
}

void SFA::printSuffix()
{
    for (SFAIDXTYPE i = 0; i < size; i++ ) 
        std::cout << i << "\t" << &concat[SA[i]] << "\n";
}


void SFA::readSFA( const char *filename )
{
    size = fio::read<SFAIDXTYPE>( SA, filename, 0 );
}

void SFA::readDoc( const char *filename )
{
    size_t n = fio::read<SFACHARTYPE>( concat, filename, 1 );
    assert( (SFAIDXTYPE)n == size );
    concat[size] = '\0';
}

void SFA::readLCP( const char *filename )
{
    int n = fio::read<LCPTYPE>( LCP, filename, 0 );
    assert( (SFAIDXTYPE)n == size );
}

void SFA::readMLCP( const char *filename )
{
    ia.load(filename, size-1);
    mLCP = ia.getIntervals();
}

void SFA::load( const char *sfa_file,
                const char *doc_file,
                const char *lcp_file,
                const char *mcp_file )
{
    readSFA( sfa_file );
    readDoc( doc_file );
    readLCP( lcp_file );
    readMLCP( mcp_file );
}

int SFA::lcp(const SFACHARTYPE* a, const SFACHARTYPE* b) 
{
  int i;
	for (i = 0; *a && *b && *a == *b; i++, a++, b++) {}
	return i;
}

void SFA::search( std::vector<SFAIDXTYPE> &pos,
                  const SFACHARTYPE *srch, 
                  int len )
{
    if ( srch[0] < concat[SA[0]] ) return;
    if ( srch[0] > concat[SA[size-1]] ) return;
        
    SFAIDXTYPE L = 0, R = size-1;
    SFAIDXTYPE M;
    while( R-L>1) {
        M = (L+R)/2;
        int cmp = strncmp( (char *)srch, (char *)&concat[SA[M]], len );
        if ( cmp <= 0 ) {
            R = M;
            if ( cmp == 0 ) pos.push_back(SA[M]);
        } else  { 
            L = M;
        }
    }
}

// No match: if left > right
// #  Match: right-left+1
BOUNDTYPE SFA::search( const SFACHARTYPE *srch, int len )
{
    return ( LCP != NULL && mLCP != NULL ) ?
        searchWithLCPs( srch, len ) :
        searchOnSFA( srch, len );
}

BOUNDTYPE SFA::searchWithLCPs( const SFACHARTYPE *srch, int len )
{
    SFAIDXTYPE left  = getLeftBoundWithLCPs( srch, len );
    SFAIDXTYPE right = getRightBoundWithLCPs( srch, len );

    return BOUNDTYPE(left-1, right-1);
}

BOUNDTYPE SFA::searchOnSFA( const SFACHARTYPE *srch, int len )
{
    SFAIDXTYPE left  = getLeftBoundOnSFA( srch, len );
    SFAIDXTYPE right = getRightBoundOnSFA( srch, len );

    return BOUNDTYPE(left-1, right-1);
}



SFAIDXTYPE SFA::getLeftBoundOnSFA(const SFACHARTYPE* pat, int len)
{
	int lLcp = lcp( concat+SA[0], pat );
	int rLcp = lcp( concat+SA[size-1], pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp >= len || pat[lLcp] <= concat[SA[0]+lLcp]) return 1;
	if (rLcp <  len && pat[rLcp] >= concat[SA[size-1]+rLcp]) return size+1;

	SFAIDXTYPE left = 1;
	SFAIDXTYPE right = size;
	while (right-left > 1)
	{
    SFAIDXTYPE mid = left + floor((right-left)/2.0);
		LCPTYPE mLcp = lLcp <= rLcp? lLcp : rLcp;
		mLcp += lcp(concat+SA[mid-1]+mLcp, pat+mLcp);
		if (mLcp >= len || pat[mLcp] <= concat[SA[mid-1]+mLcp]) { 
            right = mid; rLcp = mLcp; 
        }
		else { 
            left  = mid; lLcp = mLcp; 
        }
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

SFAIDXTYPE SFA::getRightBoundOnSFA(const SFACHARTYPE* pat, int len) 
{
	int lLcp = lcp(concat+SA[0], pat);
	int rLcp = lcp(concat+SA[size-1], pat);

	if (lLcp <  len && pat[lLcp] <= concat[SA[0]+lLcp]) return 0;
	if (rLcp >= len || pat[rLcp] >= concat[SA[size-1]+rLcp]) return size;

	SFAIDXTYPE left = 1;
	SFAIDXTYPE right = size;
	while (right-left > 1)
	{
		SFAIDXTYPE mid = left + floor((right-left)/2.0);
		LCPTYPE mLcp = lLcp <= rLcp? lLcp : rLcp;
		mLcp += lcp(concat+SA[mid-1]+mLcp, pat+mLcp);
		if (mLcp >= len || pat[mLcp] >= concat[SA[mid-1]+mLcp]) { 
            left  = mid; lLcp = mLcp; 
        }
		else { 
            right = mid; rLcp = mLcp; 
        }
	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

SFAIDXTYPE SFA::getLeftBoundWithLCPs(const SFACHARTYPE* pat, int len)
{
    if ( verbose ) std::cout << "Left bound:\n";

	int lLcp = lcp( concat+SA[0], pat );
	int rLcp = lcp( concat+SA[size-1], pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp >= len || pat[lLcp] <= concat[SA[0]+lLcp]) return 1;
	if (rLcp <  len && pat[rLcp] >= concat[SA[size-1]+rLcp]) return size+1;

	SFAIDXTYPE left  = 1;
	SFAIDXTYPE right = size;
  LCPTYPE mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        SFAIDXTYPE mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);

        if ( lLcp >= rLcp ) {
            lMcp = getMcp(left, mid);
            if ( lMcp >= lLcp )
                mLcp = lLcp + lcp(concat+SA[mid-1]+lLcp, pat+lLcp);
            else {
                mLcp = lMcp;
            }
        } else {
            rMcp = getMcp(mid,right);
            if ( rMcp >= rLcp ) 
                mLcp = rLcp + lcp(concat+SA[mid-1]+rLcp, pat+rLcp);
            else {
                mLcp = rMcp;
            }
        }

		if (mLcp >= len || pat[mLcp] >= concat[SA[mid-1]+mLcp]) { 
            right = mid; rLcp = mLcp; 
        }
		else { 
            left  = mid; lLcp = mLcp; 
        }
        if (verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);
	}

    if ( verbose ) std::cout << "left bound:" << right << "\n";
	return right;
}

SFAIDXTYPE SFA::getRightBoundWithLCPs(const SFACHARTYPE* pat, int len) 
{
    if ( verbose ) std::cout << "Right bound:\n";

	int lLcp = lcp( concat+SA[0], pat );
	int rLcp = lcp( concat+SA[size-1], pat );
    if ( verbose ) printf("lLcp:%d\trLcp:%d\n", lLcp, rLcp);

	if (lLcp <  len && pat[lLcp] <= concat[SA[0]+lLcp]) return 0;
	if (rLcp >= len || pat[rLcp] >= concat[SA[size-1]+rLcp]) return size;

	SFAIDXTYPE left = 1;
	SFAIDXTYPE right = size;
  LCPTYPE mLcp, lMcp, rMcp;
	while (right-left > 1)
	{
        lMcp = rMcp = 0;
        SFAIDXTYPE mid = (left+right)/2;

        if ( verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d lMcp:%d rMcp:%d\n", left, right, mid, lLcp, rLcp, lMcp, rMcp);
        if ( rLcp >= lLcp ) {
            rMcp = getMcp(mid,right);
            if ( rMcp >= rLcp ) 
                mLcp = rLcp + lcp(concat+SA[mid-1]+rLcp, pat+rLcp);
            else {
                mLcp = rMcp;
            }
        }
        else {
            lMcp = getMcp(left, mid);
            if ( lMcp >= lLcp ) 
                mLcp = lLcp + lcp(concat+SA[mid-1]+lLcp, pat+lLcp);
            else {
                mLcp = lMcp;
            }
        } 
		if (mLcp >= len || pat[mLcp] >= concat[SA[mid-1]+mLcp]) { 
            left  = mid; lLcp = mLcp; 
        }
		else { 
            right = mid; rLcp = mLcp; 
        }
        if ( verbose ) 
            printf("left:%lu right:%lu mid:%lu lLcp:%d rLcp:%d mLcp:%d\n", left, right, mid, lLcp, rLcp, mLcp);

	}

    if ( verbose ) std::cout << "right bound:" << left << "\n";
	return left;
}

LCPTYPE SFA::getMcp(SFAIDXTYPE l, SFAIDXTYPE r)
{
    assert( r >= l );
    if ( l == r ) return 0;
    if ( r-l == 1 ) return (LCPTYPE)LCP[r-1];
    
    SFAIDXTYPE mid = (l+r)/2;
    assert(mid<size);
    return (LCPTYPE)mLCP[mid-1];
}


char* SFA::getSuffix(SFAIDXTYPE p)
{
    assert( p <= size );
    return (char*)&concat[SA[p]];
}

SFAIDXTYPE SFA::getAt(SFAIDXTYPE p)
{
    assert( p <= size );
    return SA[p];
}






