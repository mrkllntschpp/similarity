# Similarity Metrics for Research

In today's data-driven research landscape, quantifying the similarity or difference between datasets is essential for various applications. Datasets come in diverse forms, including scalar values, 1-dimensional vectors, 2-dimensional matrices, and n-dimensional tensors. To enable this comparison, various mathematical measures are employed to assess the similarity or dissimilarity (distance) between these data structures.

## Overview

This repository provides a collection of similarity and distance metrics designed to address a wide range of data types. These metrics help researchers, data scientists, and engineers analyze, compare, and interpret their data effectively. Whether you are working with numerical data, text, images, or any other type of information, these metrics offer versatile solutions for your research needs.



### Distance and Similarity Families

The similarity metrics in this repository are organized into different families, each catering to specific use cases. Below, we provide a brief overview of these families:

- **$L_p$ Minkowski Family:** Contains measures related to the generalized formula, where p encompasses a wide range of distance norms from the $L_1$ (city block distance) to the Chebyshev $L_∞$ distance.

- **$L_1$ Family:** Contains measures related to the absolute difference $L_1$ distance (i.e., the sum of absolute differences between corresponding data points).

- **Intersection Family:** Contains measures related to the intersection of X and Y. The $\min{(X, Y)}$ or $\max{(X, Y)}$ term appears in either the denominator or numerator for this family.

- **Inner Product Family:** Contains measures related to the summed inner product or dot product of X and Y (i.e., the sum of the product of corresponding data points).

- **Fidelity (Squared-Chord) Family:** Contains measures related to the sum of the square root of the inner product, referred to as the Fidelity similarity.

- **Squared $L_2$ ($\chi^2$) Family:** Contains measures related to the square of the $L_2$ (Euclidean) norm, with some measures providing asymmetric responses if X and Y are swapped.

- **Shannon's Entropy Family:** Contains measures related to Shannon's concept of probabilistic uncertainty or entropy.

- **Combination Family:** Contains measures that combine concepts from multiple families, offering a comprehensive approach to similarity.

- **Vicissitude Family:** Contains a variety of measures.

Please note that while these metrics offer powerful ways to quantify similarity, there are some caveats to their implementation, most notably errors associated with division by zero or taking the log of zero.

## Example of Distance Metrics in Lp Minkowski Family

Here are just a few of the common distance metrics provided by this repository:

1. **$L_1$-Norm (Manhattan Distance):**
   - Formula: $\sum{|X_i - Y_i|}$
   - Description: Computes the sum of the absolute differences between corresponding data points. It is related to the mean absolute error.

2. **$L_2$-Norm (Euclidean Distance):**
   - Formula: $\sqrt{\sum{(X_i - Y_i)^2}}$
   - Description: Calculates the square root of the sum of squared differences between corresponding data points. This is related to the root mean squared error.

3. **$L_p$-Norm (Minkowski Distance):**
   - Formula: $\left({\sum{(X_i - Y_i)^p}}\right)^\frac{1}{p}$
   - Description: Calculates the difference raised to the $p$ power and then the root of $1/p$ is applied.

4. **$L_∞$-Norm (Chebyshev Distance):**
   - Formula: $\max(|X_i - Y_i|)$
   - Description: Measures the maximum distance between the two datasets.

## Getting Started

To get started with using these similarity metrics, please refer to the documentation and code examples provided in the relevant folders and files.

## Contributions

Contributions and improvements to this repository are highly encouraged. If you have a new similarity metric to add or you've found a bug in an existing one, please consider contributing by creating a pull request.

## License

This project is licensed under the [MIT License](LICENSE) - see the [LICENSE](LICENSE) file for details.
