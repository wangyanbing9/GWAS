import random


class Ant:

    @staticmethod
    def ant_first_snp(snps_sequence, snps_scores):
        """
        定义一个蚂蚁来依据snps_scores来选择1个snp
        """
        random.seed()
        x = random.uniform(0, sum(snps_scores))
        cumulative_probability = 0.0
        for item, item_probability in zip(snps_sequence, snps_scores):
            cumulative_probability += item_probability
            if x < cumulative_probability:
                return item

    @staticmethod
    def ant_second_snp(first_snp, dict_pair, dict_prob):
        """
        定义一个蚂蚁来依据条件概率dict_prob来选择第2个snp
        """
        random.seed()
        x = random.uniform(0, sum(dict_prob[first_snp]) - 1)
        cumulative_probability = 0.0
        for item, item_probability in zip(dict_pair[first_snp], dict_prob[first_snp]):
            cumulative_probability += item_probability
            if x < cumulative_probability:
                return item[1]
