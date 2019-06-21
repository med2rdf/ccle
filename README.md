# ccle
CCLE RDF converter

## Get Dataset  
Download data from the following web site .
* [CCLE](https://portals.broadinstitute.org/ccle)
```
- CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct
- CCLE_DepMap_18q3_maf_20180718.txt
- CCLE_sample_info_file_2012-10-18.txt
```

* [Cellosaurus](https://web.expasy.org/cellosaurus/)
```
- cellosaurus.txt
```

## Check docker-compose.yaml
Edit docker-compose.yaml and build docker container using docker-compose .
```
docker-compose.yaml
--------------------------------
#      - [directory downloading ccle data]:/work/data/raw         for mount
#      - [RDF file output directory]:/work/output                 for mount
--------------------------------
```

## Create docker container
```
docker-compose up -d

docker-compose exec ccle_cnvtr bash
```

## Create RDF files
Check script and raw files in container .
```
/work/
    ┣ data/
    ┃   ┗ raw/
    ┃       ┣ CCLE_DepMap_18q3_RNAseq_RPKM_20180718.gct
    ┃       ┣ CCLE_DepMap_18q3_maf_20180718.txt
    ┃       ┣ CCLE_sample_info_file_2012-10-18.txt
    ┃       ┗ cellosaurus.txt
    ┃
    ┣ scripts/
    ┃   ┣ __init__.py
    :   ┣ ccle_to_rdf.py
        ┣ omics_to_rdf.py           
        ┗ preprocessing.py
```

Run as follow . (It takes about 30 minites / 10 core)
```
cd /work/scripts
python3 preprocessing.py
python3 ccle_to_rdf.py
```
