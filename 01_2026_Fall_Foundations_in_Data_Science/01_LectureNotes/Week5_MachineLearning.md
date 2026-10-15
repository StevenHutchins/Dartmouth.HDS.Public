## Week 5 Lecture Notes - Introduction to Machine Learning

These are my notes from HSE 711: Foundations in Data Science. Week 5 was designed to provide an introduction to machine learning.

## Key Concepts

| # | Concept | Description | Example |
|---|------|--------|-------------|
| 01 | Random forest | Ask lots of simple yes/no checklists (decision trees) and go with the majority vote. | 3 "nurses" vote on a 72-year-old with 2 prior visits: Yes, Yes, No. 2 of 3 say Yes, so predict readmitted. |
| 02 | Support vector machine (SVM) | Draw the widest possible fence between two groups. The dots closest to the fence (support vectors) decide where it goes. | Blue dots = didn't come back, red = came back. A new patient lands on the red side, so predict readmitted. |
| 03 | K-nearest neighbors (KNN) | Predict from the K past cases most similar to the new one. | New patient: age 60, BMI 30. Nearest 3: (58, 29) came back, (63, 32) came back, (57, 33) did not. 2 of 3 came back, so predict readmitted. |
| 04 | Principal Component Analysis (PCA) | Squash many measurements into a few summary scores by finding what varies the most. It doesn't know the group labels. | Picks height because it varies a lot, but height doesn't separate who comes back, so the groups stay mixed. |
| 05 | Partial Least Squares Discriminant Analysis (PLS-DA) | Squash many measurements into a few summary scores by finding what best separates the groups. It uses the labels. | Picks the blood test value, which cleanly splits who came back from who didn't. |
| 06 | Elastic net | A regression formula that shrinks every factor's weight and drops the useless ones to zero. Each factor must improve predictions enough to "earn" its weight. | risk = 3×age + 1.5×BMI + 0.8×smoking + 0.2×shoe size becomes 1.67×age + 0.67×BMI + 0.2×smoking, and shoe size is dropped. |

## Packages

| #  | Package   | Description | Example |
|----|-----------|-------------|---------|
| 01 | RANN | Fast "nearest neighbor" search: finds the most similar rows to a given row.<br>caret uses it behind the scenes to fill in missing values with KNN<br>(preProcess method = "knnImpute"). | A patient is missing BMI. RANN finds the 5 most similar patients,<br>and their BMIs are used to fill in the blank. |
| 02 | tidyverse | A bundle of packages for cleaning, reshaping, and plotting data.<br>Includes dplyr (filter, select, mutate), ggplot2 (charts),<br>readr (read files), and tidyr (reshape). | data %>% filter(age > 65) %>% ggplot(aes(age, bmi)) + geom_point() |
| 03 | caret | One consistent toolkit for machine learning: splitting data, preprocessing,<br>training many model types (random forest, SVM, KNN, elastic net, PLS-DA)<br>with cross-validation, and measuring accuracy. | train(readmit ~ ., data = train_set, method = "rf",<br>trControl = trainControl(method = "cv", number = 5)) |
