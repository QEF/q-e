include(UseLATEX)

# By default the manuals are generated and installed in pdf format
# To generat the manuals use the following command: make pdf

# It is possible to generate user_guide
# in HTML format with the following command: make html

set(LATEX_DEFAULT_BUILD "pdf")

add_latex_document(
    user_guide.tex
    INPUTS quote.tex html.sty
    IMAGES quantum_espresso.png)

add_latex_document(
    constraints_HOWTO.tex
    FORCE_PDF)

add_latex_document(
    brillouin_zones.tex
    IMAGE_DIRS images
    FORCE_PDF)

add_latex_document(
    plumed_quick_ref.tex
    IMAGE_DIRS pic
    FORCE_PDF)

install(
    FILES
        ${PROJECT_BINARY_DIR}/Doc/user_guide.pdf 
        ${PROJECT_BINARY_DIR}/Doc/constraints_HOWTO.pdf
        ${PROJECT_BINARY_DIR}/Doc/brillouin_zones.pdf
        ${PROJECT_BINARY_DIR}/Doc/plumed_quick_ref.pdf
    DESTINATION 
        share/Doc)
