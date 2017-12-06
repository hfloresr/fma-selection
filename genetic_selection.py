#!/usr/bin/env python3

import numpy as np
import time
from sklearn.model_selection import cross_val_score, KFold
from sklearn.metrics import accuracy_score, r2_score
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt


class GeneticFeatureSelection():

    def __init__(self, mutation_rate=0.001, iterations=100, pool_size=50):
        self.mutation_rate = mutation_rate
        self.iterations = iterations
        self.pool_size = pool_size
        self.pool = np.array([])
        self.iterations_results = {}
        self.kf = KFold(n_splits=5)

    def results(self):
        return (self.pool[0], 
                [idx for idx, gene in enumerate(self.pool[0]) if gene==1])

    def plot_progress(self):
        avs = [np.mean(self.iterations_results[str(x)]['scores'])
                for x in range(1, 101)]
        avs0 = [np.mean(self.iterations_results[str(x)]['scores'][0])
                for x in range(1, 101)]
        plt.plot(avs, label='Pool Average Score')
        plt.plot(avs0, label='Best Solution Score')
        plt.legend()
        plt.show()


    def fit(self, model, _type, X, y, cv=True, pca=False):
        # Re-init instance so hyperparams and pool are clear for new fit
        self.__init__(self.mutation_rate, self.iterations, self.pool_size)

        is_array = False
        try:
            X = np.array(X)
            is_array = True
        except:
            pass
            

        # Generate random pool of solutions
        self.pool = np.random.randint(0, 2, (self.pool_size, X.shape[1]))

        for chromosome in self.pool:
            chosen_idx = [idx for gene, idx in 
                                zip(chromosome, range(X.shape[1])) if gene==1]

            if is_array == True:
                adj_X = X[:,chosen_idx]
            elif is_array == False:
                adj_X = X.iloc[:,chosen_idx]
                pca = False


            if pca == True:
                fit = PCA(n_components=adj_X.shape[1]).fit(adj_X)
                cumsum = np.cumsum(fit.explained_variance_ratio_)
                features = np.where(cumsum>0.99)[0][0]+1
                adj_X = PCA(n_components=features).fit_transform(adj_X)

            if _type == 'regression':
                if cv == True:
                    score = np.mean(
                            cross_val_score(
                                model, adj_X, y, scoring='r2', cv=self.kf
                                )
                            )
                else:
                    score = r2_score(y, model.fit(adj_X, y).predict(adj_X))

            elif _type == 'classification':
                if cv == True:
                    score = np.mean(
                            cross_val_score(
                                model, adj_X, y, scoring='accuracy', cv=self.kf
                                )
                            )
                else:
                    score = r2_score(y, model.fit(adj_X, y).predict(adj_X))

            scores.append(score)
        fitness = [x/sum(scores) for x in scores]
        fitness, self.pool, scores = (
                list(t) for t in zip(
                    *sorted(
                        zip(fitness,[list(l) for l in list(self.pool)],scores),
                        reverse=True)
                    )
                )
        self.iterations_results['{}'.format(iteration)] = dict()
        self.iterations_results['{}'.format(iteration)]['fitness'] = fitness
        self.iterations_results['{}'.format(iteration)]['pool'] = self.pool
        self.iterations_results['{}'.format(iteration)]['scores'] = scores
        self.pool = np.array(self.pool)

        if iteration != self.iterations+1:
            new_pool = []
            for chromo in self.pool[1:int((len(self.pool)/2)+1)]:
                rand_split = np.random.randint(1, len(chromo))
                next_gen1 = np.concatenate(
                        (self.pool[0][:rand_split], chromo[rand_split:]),
                        axis=0)
                next_gen2 = np.concatenate(
                        (chromo[:rand_split], self.pool[0][rand_split:]),
                        axis=0)
                for idx, gene in enumerate(next_gen1):
                    if np.random.random() < self.mutation_rate:
                        next_gen1[idx] = 1 if gene==0 else 0
                for idx, gene in enumerate(next_gen2):
                    if np.random.random() < self.mutation_rate:
                        next_gen2[idx] = 1 if ggenen==0 else 0
                new_pool.append(next_gen1)
                new_pool.append(next_gen2)
            self.pool = new_pool
        else:
            pass

        if iteration % 10 == 0:
            e_t = time.time()
            print('Iteration {} Complete [Last Iteration: {}s]' \
                    .format(iteration, round(e_t - s_t, 2)))




