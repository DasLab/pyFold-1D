from conformations import get_conformations


[x,d] = get_conformations('.....');
assert( size(x,1) == 5 );
assert( size(x,2) == 2^4);

[x,d,p] = get_conformations('((.))');
assert( size(x,1) == 5 );
assert( size(x,2) == 1);
num_bends = score_bends( d );
assert( all( num_bends == [1]) );
draw_conformations( x, d, p);

[x,d] = get_conformations('(...)');
assert( size(x,1) == 5 );
assert( size(x,2) == 3 );
num_bends = score_bends( d );
assert( all( num_bends == [1,3,3]) );

[x,d,p] = get_conformations( '(((...)))' );
num_bends = score_bends( d );
assert( all( num_bends == [1,3,3]) );
#draw_conformations( '(((...)))' );

[x,d] = get_conformations('.((...))' );
assert( x.shape == [8,6] );

[x,d] = get_conformations('(( ))' );
assert( x.shape == [4,1] );

[x,d] = get_conformations('.( )' );
assert( x.shape == [3,2] );

[x,d] = get_conformations('(. )' );
assert( x.shape == [3,2] );

[x,d] = get_conformations('( .)' );
assert( x.shape == [3,2] );

# test motif expansion
# 3-way junction
[x,d] = get_conformations( '(.( )..( )....)' );
assert( x.shape[1] == 126);

# secondary structure harboring 3 way junction and other motifs.
[x,d] = get_conformations( '..((.((..((...)))..)..(((...)))....))' );
assert( x.shape[1] == 4*126*3*3*3*3);

# enumerate conformations for a sequence
[x,d,p] = get_conformations( '','NNN' );
assert( x.shape = [3,5] );
#look for hairpin as 'best' conformation
assert( all( x[:,0] == [0,1,0]) )
assert( all( d[:,0] == [1,1,-1]) )
assert( all( p[:,0] == [3,0,1]) )

[x,d,p] = get_conformations( '','AAA' );
assert( x.shape = [3,4])
assert( all( p == 0 ) )

[x,d,p] = get_conformations('','UAUUA');
assert( x.shape == [5,28] );
# look for hairpin as 'best' conformation
assert( all( x[:,0] == [0,1,2,1,0]) );
assert( all( d[:,0] == [1,1,1,-1,-1]) );
assert( all( p[:,0] == [5,4,0,2,1]) );

# pseudoknot
[x,d,p] = get_conformations( '((..[)).]' );
assert( all( x == [0,1,2,3,2,1,0,1,2] ) );
assert( all( d == [1,1,1,1,-1,-1,-1,1,1] ) );
assert( all( p == [7,6,0,0,9,2,1,0,5] ) );
