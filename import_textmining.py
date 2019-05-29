import textmining
from pymongo import MongoClient


def main():
    output = textmining.main([["anti-glycenia", "blood sugar", "Daf16", "insulin", "foxo", "survival", "lifespan",
                               "alc"],
                              ["obesity", "cholesterol", "HDL", "LDL", "lipoprotein", "liver", "pancreas",
                               "islands of langerhands", "p cells"],
                              ["gallic acid", "charetin", "vicine", "polypeptide p", "charetin"]])

    client = MongoClient("mongodb+srv://textminingdata-m49re.azure.mongodb.net/test", 27017,
                         username="admin", password="blaat1234")

    db = client.gourd_goru
    collection = db.miningResults

    # collection.insert_many(output)


main()
