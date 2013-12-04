\name{misclassify}
\alias{misclassify}
\title{Adds some misclassification error to binary data}
\usage{
  misclassify(binary.vector, error.1.0 = 0.05,
    error.0.1 = 0.05)
}
\arguments{
  \item{binary.vector}{a vector of binary values.}

  \item{error.1.0}{1 to 0 misclassification rate.}

  \item{error.0.1}{0 to 1 misclassification rate.}
}
\value{
  a binary vector
}
\description{
  Introduces some misclassification in binary data. The
  number of values misclassified is determined by the
  misclassification rates.
}
\examples{
{

# Example 1: a vector of binaries coded '0' and '1'
# simulate a vector of 0s and 1s
v1 <- rbinom(100,1,0.4)
# 1 to 0 misclassification rate
error.1.0 <- 0.2
# 0 to 1 misclassification rate
error.0.1 <- 0.2
# randomly misclassify the vector v1
v1.new <- misclassify(v1, error.1.0, error.0.1)

# Example 2: a vector of binaries coded '2' and '3'
# simulate a vector of 2s and 3s
v2 <- ifelse(runif(100, 0, 1) <= 0.4, 2,3)
# 1 to 0 misclassification rate
error.1.0 <- 0.2
# 0 to 1 misclassification rate
error.0.1 <- 0.2
# randomly misclassify the vector v
v2.new <- misclassify(v2, error.1.0, error.0.1)
}
}
\author{
  Amadou Gaye
}

