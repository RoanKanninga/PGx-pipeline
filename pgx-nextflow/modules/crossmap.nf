process crossmap {

!!!!!	modules Crossmap

    tag 'crossmap'
    input:
    tuple path(bed), path(bim), path(fam)
    path chain_file
    
    output:
    tuple path('crossmapped_plink.bed'), path('crossmapped_plink.bim'), path('crossmapped_plink.fam')
    template 'crossmap.sh'
	}