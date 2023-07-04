
import logging
import numpy as np


logging.basicConfig(filename=os.path.join(args.output_dir, 
                                          "predictions.log"),
                    format="%(asctime)s %(levelname)s : %(message)s",
                    level = logging.INFO)
def run():
    logging.INFO("Start: effect size inference")
    logging.INFO(f"vcf file:{vcf}")
    logging.INFO(f"expression file:{par_file}")

    print(f"Successfully found {__file__}.")

    logging.INFO(f"output written to:{output_dir/predictions.bed}")
    logging.INFO("Finished")

    raise NotImplementedError
