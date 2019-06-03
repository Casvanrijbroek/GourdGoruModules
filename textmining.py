from nltk.corpus import stopwords
from pymongo import MongoClient
from itertools import product
from Bio import Entrez
import nltk
import sys
import re

""" This script can be used to perform textmining on PubMed articles based on a set of search terms. The search terms 
need to be given in multiple categories (nested lists), because the script only compares keywords of different 
categories with each other. To use this script, simply import it into your application and pass the nested list of
search terms to the main function (other functions can be called independently to fit a different pipeline too). 

Author: Cas van Rijbroek
Date: 03-jun-2019
Version: 1.2
"""


def textmine(search_words):
    """ This function can be called from another script to perform textmining on search words.

    :param search_words: A nested list of search terms, each list should contain words that are related to each
    other and therefor need not be compared (such as terms related to a single disease).
    :return: A list of dictionaries that contain the found co-occurrence between the given search terms. For each
    article found some general information is stored together with the assigned co-occurrence score.
    """
    products = create_products(search_words)

    return main(products)


def textmine_new_term(new_term, upper_term):
    """ This function is called when the script is called from the command line. The intent behind this function is to
    add new terms to the existing pool of results in the database.

    :param new_term: A string containing a new search word.
    :param upper_term: A string containing the upper term of the search word or "new" if the word is part of a new
    category.
    :return: A list of dictionaries that contain the found co-occurrence between the given search terms. For each
    article found some general information is stored together with the assigned co-occurrence score.
    """
    search_words = get_search_words(new_term, upper_term)
    products = create_products_additional_term(search_words)

    return main(products)


def main(products):
    """ This is the main function that calls the other functions required to perform textmining on products.

    :param products: A list with sets that contain every wanted combination of search words.
    :return: A list of dictionaries that contain the found co-occurrence between the given search terms. For each
    article found some general information is stored together with the assigned co-occurrence score.
    """
    Entrez.email = "casvanrijbroek@hotmail.com"
    queries = create_queries(products)
    print(queries)
    output = mine(list(products), queries)
    print(output)

    return output


def get_search_words(new_term, upper_term):
    """ This function fetches all existing search terms from the database with the exclusion of the ones contained
    in the same category as the new term. The new term is also included in the output.

    :param new_term: A string containing a new search word.
    :param upper_term: A string containing the upper term of the search word or "new" if the word is part of a new
    category.
    :return: A nested list of search terms, each list should contain words that are related to each
    other and therefor need not be compared (such as terms related to a single disease).
    """
    search_words = []

    client = MongoClient("mongodb+srv://textminingdata-m49re.azure.mongodb.net/test", 27017,
                         username="admin", password="blaat1234")
    db = client.gourd_goru
    collection = db.searchWords

    if upper_term == "new":
        terms = collection.find({})
        search_words.append([new_term])
    else:
        terms = collection.find({"upper_term": {"$ne": upper_term}})
        search_words.append([new_term])

    for category in terms:
        search_words.append(category["sub_terms"])

    return search_words


def create_products_additional_term(search_words):
    """ This function creates all wanted combinations of search words based on the addition of a new search word into
    the database. This means that the new search word will be combined with all other search words from different
    categories.

    :param search_words: A nested list of search terms, each list should contain words that are related to each
    other and therefor need not be compared (such as terms related to a single disease).
    :return: A list with sets that contain every wanted combination of search words.
    """
    products = set()

    for category in search_words:
        if category != search_words[0]:
            for word in category:
                products.add((word, search_words[0][0]))

    return products


def create_products(search_words):
    """ This function creates combinations of the search words. It combines every word with each other with the
    exception of words from the same list.

    :param search_words: A nested list of search terms, each list should contain words that are related to each
    other and therefor need not be compared (such as terms related to a single disease).
    :return: A list with sets that contain every wanted combination of search words.
    """
    products = set()

    for i in range(len(search_words) - 1):
        for i2 in range(i, len(search_words) - 1):
            products.update(set(product(search_words[i], search_words[i2 + 1])))

    return products


def create_queries(products):
    """ Transforms a list of search word combinations into queries suitable for PubMed searches.

    :param products: A list with sets that contain every wanted combination of search words.
    :return: A list containing PubMed queries.
    """
    queries = []

    for two in products:
        queries.append(" AND ".join(two))

    return queries


def mine(search_words, queries):
    """ This is the main text mining function. It calls other functions in this script to perform PubMed searches using
    queries and assigns co-occurrence scores to the returned articles.

    :param search_words: A nested list that contains every wanted combination of search words.
    :param queries: A list containing PubMed queries (this list must be complimentary to the list of search words).
    :return: A list of dictionaries that contain the found co-occurrence between the given search terms. For each
    article found some general information is stored together with the assigned co-occurrence score.
    """
    full_output = []

    for index in range(len(search_words)):
        ids, mesh_terms = get_pubmed_ids(queries[index])
        output = []

        if len(ids) > 0:
            articles = get_articles(ids)
            abstracts = preprocess_data(articles)

            for article_index in range(len(abstracts)):
                if abstracts[article_index] is False:
                    score = 0
                else:
                    score = calculate_score(mesh_terms, abstracts[article_index],
                                            articles["PubmedArticle"][article_index]["MedlineCitation"]
                                            ["Article"]["ArticleTitle"].lower())

                output.append(format_output(articles["PubmedArticle"][article_index], score))

        full_output.append({"keywords": search_words[index],
                            "articles": output})

    return full_output


def get_pubmed_ids(query):
    """ This function performs a PubMed search using a single query and returns the list of PubMed ID's and related
    mesh terms.

    :param query: A single PubMed query
    :return: A list of PubMed ID's returned from the Search (this list can be empty),
    Mesh terms related to the given search terms (if NCBI didn't identify mesh terms the original search word(s) will
    be located in this list instead.
    """
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
    """ This function fetches articles from PubMed based on a list of PubMed ID's.

    :param pubmed_ids: A list of PubMed ID's.
    :return: A PubMed record (dictionary) containing all found articles.
    """
    handle = Entrez.efetch(db="pubmed", id=",".join(pubmed_ids), retmode="xml")
    record = Entrez.read(handle)
    handle.close()

    return record


def preprocess_data(data):
    """ Preprocesses abstracts (or any form of text) by tokenizing it and removing unwanted words to speed up the
    process of textmining.

    :param data: A PubMed record (dictionary).
    :return: A nested list containing tokenized (lists with words (strings)) abstracts or a False boolean if no
    abstract was found for a certain article.
    """
    abstracts = []
    unwanted_words = ["PROTOCOL", "ABSTRACT", "INTRODUCTION", "MATERIAL", "METHODS", "CONCLUSION", "PATIENTS", "aim",
                      "of", "study", "although", "well", "include", "results"]

    for article in data["PubmedArticle"]:
        if "Abstract" in article["MedlineCitation"]["Article"]:
            abstract = ""

            for line in article["MedlineCitation"]["Article"]["Abstract"]["AbstractText"]:
                abstract += line

            tokens = nltk.word_tokenize(abstract)
            stop_words = stopwords.words("english")
            tokens = [token for token in tokens if token.isalpha()]
            tokens = [token.lower() for token in tokens if token.lower() not in stop_words]
            tokens = [tokens for token in tokens if token not in unwanted_words]
            abstracts.append(tokens)

        else:
            abstracts.append(False)

    return abstracts


def calculate_score(search_words, abstract, title):
    """ This function calculates the co-occurrence score of a single pre-processed abstract (or any form of text).
    Occurrence in the abstract is worth 5 points while occurrence in the title is 10 points.

    :param search_words: A nested list that contains every wanted combination of search words.
    :param abstract: A single tokenized abstract (list with words (strings)).
    :param title: The title of the article belonging to the abstract.
    :return: An integer indicating the co-occurrence score of the given abstract + title).
    """
    score = 0

    for words in search_words:
        for word in words:
            if word in title:
                score += 10

            score += abstract.count(word) * 5

    return score


def format_output(article, score):
    """ This function formats the information of a PubMed article and a co-occurrence score into a dictionary format
    suitable for storage.

    :param article: A PubMed article (dictionary).
    :param score: An integer indicating the co-occurrence score of the given abstract + title).
    :return: A dictionary containing the PubMed ID, title, authors, journal, co-occurrence score and publication date
    of the given article.
    """
    medline = article["MedlineCitation"]
    pubmed_data = article["PubmedData"]
    authors = []

    if "AuthorList" in medline["Article"]:
        for author in medline["Article"]["AuthorList"]:
            if "ForeName" in author:
                authors.append(author.get("ForeName") + " " + author.get("LastName"))
            else:
                authors.append(author.get("LastName"))
    else:
        authors.append("Authors unknown")

    pub_date = "unknown"

    for date in pubmed_data["History"]:
        if date.attributes.get("PubStatus") == "pubmed":
            pub_date = "{}-{}-{}".format(date.get("Year"), date.get("Month"), date.get("Day"))

    return ({"pubmed_id": medline["PMID"],
             "title": medline["Article"]["ArticleTitle"],
             "authors": authors,
             "journal": medline["Article"]["Journal"]["Title"],
             "score": score,
             "pub_date": pub_date})


def save_results(results):
    """ This function saves the textmining results in the database.
    (NOTE: this function is NOT automatically called by the others functions. To save results this function needs to be
    specifically called by the user after textmining. The command line pipeline DOES automatically call this function.)

    :param results: A list of dictionaries that contain the found co-occurrence between the given search terms. For each
    article found some general information is stored together with the assigned co-occurrence score.
    """
    client = MongoClient("mongodb+srv://textminingdata-m49re.azure.mongodb.net/test", 27017,
                         username="admin", password="blaat1234")

    db = client.gourd_goru
    collection = db.miningResults

    collection.insert_many(results)


if __name__ == "__main__":
    """ When the script is called from the command line. This code ensures the new term is going to be textmined and
    added to the database.
    """
    textmining_results = textmine_new_term(sys.argv[1], sys.argv[2])
    save_results(textmining_results)
