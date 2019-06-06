START=$1
END=$2

for i in $(seq $START $END); do
    WT=T7_$(date +'%m%d%y')_strong3_5_w$i;
    REC=T7_$(date +'%m%d%y')_strong3_5_w${i}_r;
    python3 ./scripts/python/parse_genbank.py -w $i > ./params/$WT.yml;
    python3 ./scripts/python/parse_genbank.py -w $i -r > ./params/$REC.yml;
    echo "pinetree_run.py ./params/$WT.yml -o ./runs/$WT";
    echo "pinetree_run.py ./params/$REC.yml -o ./runs/$REC";
done;
