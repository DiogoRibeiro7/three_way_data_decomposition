simula2 <- function(N, G, nrep, dgp, ns = 100) {
    # Constants and setup
    J <- 20
    Q <- 5
    K <- 5
    R <- 2
    eps <- 1e-6

    # Result storage
    ari <- matrix(0, ns, 3)
    Xt <- array(0, c(N, J * K, ns))
    Utruet <- array(0, c(N, G, ns))

    # Simulate datasets
    for (sa in 1:ns) {
        set.seed(sa + 13)
        pg <- runif(G)
        pg <- pg / sum(pg)
        Sv <- matrix(rnorm(J * J), ncol = J)
        Sv <- t(Sv) %*% Sv
        Sv[(Q + 1):J, 1:Q] <- 0
        Sv[1:Q, (Q + 1):J] <- 0

        So <- matrix(rnorm(K * K), ncol = K)
        So <- t(So) %*% So
        So[(R + 1):K, 1:R] <- 0
        So[1:R, (R + 1):K] <- 0

        S <- kronecker(So, Sv)
        Sf <- 0.6 * matrix(rnorm(J * K * J * K), J * K)
        Sf <- t(Sf) %*% Sf
        Sf <- (S != 0) * Sf

        B <- diag(J)[, 1:Q]
        C <- diag(K)[, 1:R]
        Eta <- matrix(rnorm(Q * R * G), ncol = G)
        Mu <- t(20 * kronecker(C, B) %*% Eta)
        Muf <- 20 * (Mu != 0) * matrix(rnorm(G * J * K), ncol = J * K)

        if (dgp == 2) Mu <- Muf
        if (dgp == 3) S <- Sf
        if (dgp == 4) {
            Mu <- Muf
            S <- Sf
        }

        # Assuming genmixhet is defined
        list(X, Utrue, z) <- genmixhet(N, pg, Mu, replicate(G, S, simplify = FALSE))
        X <- preproa(X, 1) # Assuming preproa is a preprocessing function

        # Model fitting loop
        for (model in 1:3) {
            blike <- -Inf
            for (rep in 1:nrep) {
                set.seed(10 * sa + rep)
                U <- matrix(runif(N * G), nrow = N)
                U <- sweep(U, 1, rowSums(U), "/")

                if (model == 1) { # Assuming model 1 is kmeans
                    kmeans_result <- kmeans(X, centers = G, nstart = 10)
                    Ur <- kmeans_result$cluster
                    like <- kmeans_result$tot.withinss # Use total within-cluster sum of squares as a 'likelihood'
                }
                if (model == 2) { # Assuming model 2 is T3MIXS
                    # Placeholder implementation: initializing random clusters
                    Ur <- sample(1:G, size = N, replace = TRUE)
                    TB <- matrix(runif(J * Q), ncol = Q)
                    TC <- matrix(runif(K * R), ncol = R)
                    # A dummy 'like' value, typically you would have a real model fitting here
                    like <- -runif(1, 100, 200) # Dummy 'likelihood' for placeholder
                }
                if (model == 3) { # Assuming model 3 is T2MIXT
                    # Placeholder implementation: initializing random clusters
                    Ur <- sample(1:G, size = N, replace = TRUE)
                    TB <- orth(matrix(runif(J * K * Q * R), ncol = Q * R)) # Orthogonal basis
                    # A dummy 'like' value, again, you would replace this with actual fitting
                    like <- -runif(1, 100, 200) # Dummy 'likelihood' for placeholder
                }

                # Assuming mrand and ftoh are defined
                if (blike < like) {
                    blike <- like
                    bU <- U
                }
            }
            ari[sa, model] <- mrand(ftoh(Utrue) %*% t(ftoh(bU)))
            Xt[, , sa] <- X
            Utruet[, , sa] <- Utrue
        }
    }

    # Return results
    list(ari = ari, Xt = Xt, Utruet = Utruet)
}
