from Bio import Entrez
from itertools import product
from nltk.corpus import stopwords
import nltk
import sys
import re

Entrez.email = "casvanrijbroek@hotmail.com"


def main(search_words):
    products = create_products(search_words)
    queries = create_queries(products)
    mine(list(products), queries)


def mine(search_words, queries):
    full_output = {}

    for index in range(len(search_words)):
        ids, mesh_terms = get_pubmed_ids(queries[index])
        output = {}

        if len(ids) > 0:
            articles = get_articles(ids)
            abstracts = preprocess_data(articles)

            for article_index in range(len(abstracts)):
                if abstracts[article_index] is False:
                    score = 0
                else:
                    score = calculate_score(mesh_terms, abstracts[article_index])

                article = articles["PubmedArticle"][article_index]["MedlineCitation"]
                authors = []

                for author in article["Article"]["AuthorList"]:
                    if "ForeName" in author:
                        authors.append(author.get("ForeName") + " " + author.get("LastName"))
                    else:
                        authors.append(author.get("LastName"))

                pub_date = "unknown"

                for date in articles["PubmedArticle"][article_index]["PubmedData"]["History"]:
                    if date.attributes.get("PubStatus") == "pubmed":
                        pub_date = "{}-{}-{}".format(date.get("Year"), date.get("Month"), date.get("Day"))

                output.update({"pubmed_id": article["PMID"],
                               "title": article["Article"]["ArticleTitle"],
                               "authors": authors,
                               "journal": article["Article"]["Journal"]["Title"],
                               "score": score,
                               "pub_date": {pub_date}
                               })


def create_products(search_words):
    products = set()

    for i in range(len(search_words) - 1):
        for i2 in range(i, len(search_words) - 1):
            products.update(set(product(search_words[i], search_words[i2 + 1])))

    return products


def create_queries(products):
    queries = []

    for two in products:
        queries.append(" AND ".join(two))

    return queries


def get_pubmed_ids(query):
    mesh_terms = []

    handle = Entrez.esearch(db="pubmed", term=query, retmode="xml", retmax=100, sort="relevance")
    record = Entrez.read(handle)
    handle.close()

    for mesh in record["TranslationSet"]:
        mesh_terms.append(set([w.lower().strip('""') for w in re.findall(r'"\w+"', mesh["To"])]))

    if "ErrorList" in record:
        if "PhraseNotFound" in record["ErrorList"]:
            for phrase in record["ErrorList"]["PhraseNotFound"]:
                mesh_terms.append([phrase.lower()])

    return record["IdList"], mesh_terms


def get_articles(pubmed_ids):
    handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    return record


def preprocess_data(data):
    abstracts = []

    for article in data["PubmedArticle"]:
        if "Abstract" in article["MedlineCitation"]["Article"]:
            abstract = ""

            for line in article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]:
                abstract += line

            tokens = nltk.word_tokenize(abstract)
            stop_words = stopwords.words("english")
            tokens = [token for token in tokens if token.isalpha()]
            tokens = [token.lower() for token in tokens if token.lower() not in stop_words]
            abstracts.append(tokens)

        else:
            abstracts.append(False)

    return abstracts


def calculate_score(search_words, abstract):
    score = 0

    for words in search_words:
        for word in words:
            score += abstract.count(word) * 5

    return score


# sys.argv[0]
main([["anti-glycenia", "blood sugar", "Daf16"],
      ["gallic acid", "charetin", "vicine", "polypeptide p"],
      ["cushing's syndrome", "Diabetic Acidosis"]])
