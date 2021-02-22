#Case 4
EXOMISER=exomiser-cli-12.1.0.jar
EXOMISER_DIR=/mnt/common/WASSERMAN_SOFTWARE/Exomiser/exomiser-cli-12.1.0/
CONFIG_YML=Case4_Exomiser.yml

java -Xmx4g -jar -Djava.io.tmpdir=$PWD $EXOMISER_DIR/$EXOMISER \
         --analysis $CONFIG_YML \
        --spring.config.location=$EXOMISER_DIR
