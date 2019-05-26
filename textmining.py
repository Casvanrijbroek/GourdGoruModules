from Bio import Entrez
from itertools import product
from nltk.corpus import stopwords
import nltk
import sys

Entrez.email = "casvanrijbroek@hotmail.com"


def main(search_words):
    products = create_products(search_words)
    queries = create_queries(products)
    article_ids = get_pubmed_ids(queries)
    articles = get_articles(article_ids)
    abstracts = preprocess_data(articles)
    print(abstracts)


def create_products(search_words):
    products = set()

    for i in range(len(search_words)-1):
        for i2 in range(i, len(search_words)-1):
            products.update(set(product(search_words[i], search_words[i2+1])))

    return products


def create_queries(products):
    queries = []

    for two in products:
        queries.append(" AND ".join(two))

    return queries


def get_pubmed_ids(queries):
    records = []

    for query in queries:
        handle = Entrez.esearch(db="pubmed", term=query, retmode="xml", retmax=100, sort="relevance")
        record = Entrez.read(handle)
        records.append(record)
        handle.close()

    return records


def get_articles(pubmed_ids):
    records = []
    for search in pubmed_ids:
        if len(search["IdList"]) > 0:
            handle = Entrez.efetch(db="pubmed", id=",".join(search["IdList"]), retmode="xml")
            record = Entrez.read(handle)
            records.append(record)
            handle.close()
        else:
            records.append(False)

    return records


def preprocess_data(data):
    output = []
    index = 0

    for search in data:
        output.append([])
        abstracts = []

        if search:
            for article in search["PubmedArticle"]:
                abstract = ""
                if "Abstract" in article["MedlineCitation"]["Article"]:
                    for line in article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]:
                        abstract += line

                    abstracts.append(abstract)

            if len(abstracts) > 0:
                for abstract in abstracts:
                    tokens = nltk.word_tokenize(abstract)
                    tokens = [token for token in tokens if token.isalpha()]
                    stop_words = stopwords.words("english")
                    tokens = [token.lower() for token in tokens if token.lower() not in stop_words]
                    output[index].append(tokens)

        index += 1

    return output


# sys.argv[0]
main([["anti-glycenia", "blood sugar", "Daf16"],
      ["gallic acid", "charetin", "vicine", "polypeptide p"],
      ["cushing's syndrome", "Diabetic Acidosis"]])
