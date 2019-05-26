# SURF-2018

An effort from a little volunteer of SURF-2018 project. -- Aug 2018

## What is SURF?

The Summer Undergraduate Research Fellowships (SURF) programme offers undergraduate students the chance to work on research projects in collaboration with academic supervisors, during the summer holidays in XJTLU.
For more information, see [here](https://www.xjtlu.edu.cn/en/research/summer-undergraduate-research-fellowships).

## Research Title

Predict the target specificity of RNA methylation reader, writer and eraser (a machine learning/classification project)

## Research Abstract

RNA methylation has emerged as an important layer of gene regulation. Tens of thousands of RNA methylation sites are regulated by a small number of its regulators, including the readers (METTL3 and METTL16), erasers (FTO and ALKBH5) and readers (YTHDF1, YTHDF2 and YTHDF3). The proposed project seek to computational predict the substrate specificity of readers, writers and erasers so as to understand the regulatory mechanism of the RNA methylome. We will try to use different machine learning approaches to make the prediction and validate the results with cross-validation.

## My Effort:

This is my first time to get to machine learning, R programming, and even bioinformatics. I have learned loads of new knowledge (but I know they are just a little for a real bioinformatician) in this field within a quite short period, though the result is not satisfied, I know both success and failure are what I want, it is that taking effort as I may matters so much for me. 

- **Learn basic concept of machine learning.**

There are some sources I refered to:

1. [Machine learning](https://en.wikipedia.org/wiki/Machine_learning)
2. [Cross validation](https://en.wikipedia.org/wiki/Cross-validation_(statistics))
3. [Gradient boosting machine (gbm)](https://en.wikipedia.org/wiki/Gradient_boosting)
4. [Random forest (rf)](https://en.wikipedia.org/wiki/Random_forest)
5. [Support-vector machine (svm)](https://en.wikipedia.org/wiki/Support-vector_machine)

- **Learn application of machine learning method in bioinformatics (focusing on epigenetics such as posttranscriptional modification of RNAs).**

I detailedly read some papers that are closely related to this project, you can find them in "**Paper**" file:

1. Dynamic transcriptomic m6A decoration writers, erasers, readers and functions in RNA metabolism
2. MethyRNA A web-server for identification of N6-methyladenosine sites
3. SRAMP prediction of mammalian N6 methyladenosine (m6A) sites based on sequence derived features

- **Learn R programming, visualization and caret package.**

I read these two books recommended by Dr. Meng, you can find them in "**R programming books**" file:
1. An introduction to R
2. Use R

For machine learning implementation in R, Caret package is widely used. 

To learn how to use [Caret package](https://cran.r-project.org/web/packages/caret/index.html), I read [The Caret package Github pages](http://topepo.github.io/caret/index.html), and I took some notes when I learning Caret package which are in "**caret note**" file, as well as a paper "**Building Predictive Models in R Using the caret Package**" and a **concise instruction of caret functions**. 

Here are some other online sources that helped me quickly pick up R programming:(Of course there are so many sources online for you to explore and learn by yourself!)

1. [Tutorials for learning R | R-bloggers](https://www.r-bloggers.com/how-to-learn-r-2/)
2. [R language for programmers](https://www.johndcook.com/blog/r_language_for_programmers/)
3. [R Tutorial for Beginners: Learning R Programming - Guru99](https://www.guru99.com/r-tutorial.html)
4. [Online learning R-studio](https://www.rstudio.com/online-learning/#r-programming)

- **Work on this project**

My part was to predict the target specificity of RNA methylation writers —— METTL3， METTL14 and METTL16. I used rf to train the data, compared the results of cross validation by gbm, rf and svm, selected the top several important features and rebuilt the model, frustratingly, I got an unsatisfied result with pretty low AUC performance. The code is in "**code**" file.

## Reflection and Future

After reflection, I know I lack the basic understanding of this field for I was a year 1 student and just a beginner. Any real success in science need effort and time, and it is definitely impossible to get a good result for one that just learned for three months. Programming, data analysis, mathematics of machine learning and bioinformatics are what I need to spend time, maybe for years, to learn. To my delight, they are exactly what I am interested in.



