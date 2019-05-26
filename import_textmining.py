import textmining


def main():
    output = textmining.main([["anti-glycenia", "blood sugar", "Daf16"],
                              ["gallic acid", "charetin", "vicine", "polypeptide p"],
                              ["cushing's syndrome", "Diabetic Acidosis"]])

    print(output)


main()
