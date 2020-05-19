import pandas as pd
import numpy as np
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import RepeatedStratifiedKFold
import warnings
from sklearn import preprocessing
from sklearn.metrics import precision_score, recall_score, roc_auc_score, roc_curve, f1_score
from sklearn.metrics import precision_score, recall_score, roc_auc_score, roc_curve, accuracy_score
from sklearn import linear_model
from sklearn.linear_model import Ridge, RidgeCV, ElasticNet, LassoCV, LassoLarsCV
from tqdm import tqdm

warnings.filterwarnings('ignore')


#Concatenate Normal and Cancer matrices and retrieve the labels of patients
def prepare_data(N, C):
    # Put the patients on the rows and the genes on the columns
    N = N.T
    C = C.T
    # Put together normal patients and cancer ones
    full = pd.concat([N, C])
    # Change the name of the columns in numbers
    full.columns = range(len(full.columns))
    # Create the labels
    labels = pd.DataFrame(0, columns=['Normal', 'Cancer'], index=full.index)
    labels['Normal'][0:len(N)] = 1
    labels['Cancer'][len(N):(len(N) + len(C))] = 1
    full = pd.DataFrame(full)
    return full, labels

#It allows to extract a submatrix from N and C with same dimensions of N1 and C1
def same_number(N, C, N1, C1):
    full, labels = prepare_data(N, C)
    full_1, labels = prepare_data(N1, C1)
    sample = full.sample(len(full_1.T), axis=1)
    return pd.DataFrame(sample), labels


#Extraction of important variables by means of lasso
def lasso(N, C):
    full, labels = prepare_data(N, C)
    y = labels['Normal'].values
    clf = linear_model.Lasso(alpha=0.1)
    clf.fit(full, y)
    clf.predict(full)
    clf.score(full, y)

    ###Model Lasso regression
    model_lasso = LassoCV(alphas=[1, 0.1, 0.001, 0.0005]).fit(full, y)

    ###Model Lasso regression
    coef = pd.Series(model_lasso.coef_, index=full.columns)
    coef_selected = coef.iloc[coef.nonzero()]
    N = N.T
    C = C.T
    N.columns = range(len(N.columns))
    C.columns = range(len(C.columns))
    lasso_norm = N[list(coef_selected.index)]
    lasso_canc = C[list(coef_selected.index)]

    return lasso_norm.T, lasso_canc.T

#Extraction of important variables by means of lasso starting from a random set extracted from N and C of same cardinality of N1 and C1
def lasso_sample(N, C, N1, C1):
    full, labels = same_number(N, C, N1, C1)
    y = labels['Normal'].values
    clf = linear_model.Lasso(alpha=0.1)
    clf.fit(full, y)
    clf.predict(full)
    clf.score(full, y)

    ###Model Lasso regression
    model_lasso = LassoCV(alphas=[1, 0.1, 0.001, 0.0005]).fit(full, y)

    ###Model Lasso regression
    coef = pd.Series(model_lasso.coef_, index=full.columns)
    coef_selected = coef.iloc[coef.nonzero()]
    N = N.T
    C = C.T
    N.columns = range(len(N.columns))
    C.columns = range(len(C.columns))
    lasso_norm = N[list(coef_selected.index)]
    lasso_canc = C[list(coef_selected.index)]

    return lasso_norm.T, lasso_canc.T


class RF:
    def __init__(self, N, C, N_fused, C_fused, N_de, C_de, tumor):
        self.normal = N
        self.cancer = C
        self.fused_norm = N_fused
        self.fused_canc = C_fused
        self.de_norm = N_de
        self.de_canc = C_de
        self.tumor = tumor
        self.lasso_fused_norm, self.lasso_fused_canc = lasso(self.fused_norm, self.fused_canc)
        self.lasso_fused_norm.to_csv('./Extracted/lasso_norm' + str(tumor) + '.csv', sep=';')
        self.lasso_fused_canc.to_csv('./Extracted/lasso_canc' + str(tumor) + '.csv', sep=';')
        pd.DataFrame(N.index[self.lasso_fused_norm.index]).to_csv('./Extracted/genes_comm_lasso_' + str(tumor) + '.csv')
        self.lasso_de_norm, self.lasso_de_canc = lasso(self.de_norm, self.de_canc)
        self.lasso_random_values = self.RF_classifier_random(self.normal, self.cancer, self.fused_norm, self.fused_canc, self.tumor)
        self.lasso_de_values = self.RF_classifier(self.lasso_de_norm, self.lasso_de_canc, self.tumor, 'DE')
        self.extracted_values = self.RF_classifier(self.lasso_fused_norm, self.lasso_fused_canc, self.tumor,'IC')

    #One hundred times normal/cancer classification of samples
    # It classify a set of features N1 and C1
    def RF_classifier(self, N1, C1, tumor, what):
        n_reps = 100
        auc_scores = np.zeros(n_reps)
        accuracy_scores = np.zeros(n_reps)
        f1_scores = np.zeros(n_reps)
        auc_scores1 = np.zeros(n_reps)
        accuracy_scores1 = np.zeros(n_reps)
        f1_scores1 = np.zeros(n_reps)
        higher = 0
        lower = 0
        maximum = 0
        minimum = 100
        for u in range(0, n_reps):
            fused, labels = prepare_data(N1, C1)

            # prepare the data
            X = fused
            X = X.values
            y = labels
            y = y['Normal'].values

            # save all the predictions
            y_train_all = []
            y_test_all = []
            y_pred_all = []
            y_predictions_all = []
            train_probs_all = []
            train_predictions_all = []

            # cross validation
            skf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=50)
            for train_index, test_index in skf.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]

                # scale the data
                scaler = preprocessing.StandardScaler().fit(X_train)
                X_train = scaler.transform(X_train)
                X_test = scaler.transform(X_test)

                # Create the random forest instance
                mod = RandomForestClassifier(n_estimators=20, max_features=0.4, max_depth=10, min_samples_leaf=2)
                # Fit the data
                mod.fit(X_train, y_train)
                # Predict the data
                y_pred = [x[1] for x in mod.predict_proba(X_test)]
                y_predictions = mod.predict(X_test)
                train_predictions = mod.predict(X_train)
                train_probs = [x[1] for x in mod.predict_proba(X_train)]


                # Save all the predictions
                y_pred_all += list(y_pred)
                y_predictions_all += list(y_predictions)
                y_train_all += list(y_train)
                y_test_all += list(y_test)
                train_predictions_all += list(train_predictions)
                train_probs_all += list(train_probs)

            auc_scores[u] = roc_auc_score(y_test_all, y_pred_all)
            accuracy_scores[u] = accuracy_score(y_test_all, y_predictions_all)
            f1_scores[u] = f1_score(y_test_all, y_predictions_all)

        values = pd.DataFrame(columns=['auc', 'accuracy',  'f1'])
        values['auc'] = auc_scores
        values['accuracy'] = accuracy_scores
        values['f1'] = f1_scores
        values.to_csv('./Auc_acc_f1/auc_acc_f1_'+str(what)+'_' + str(tumor) + '.csv')

        return values

    #It classify a set of features extracted with lasso from a random set of genes starting from the same dimensions of N1 and C1
    def RF_classifier_random(self, N, C, N1, C1, tumor):
        n_reps = 100
        auc_scores = np.zeros(n_reps)
        accuracy_scores = np.zeros(n_reps)
        f1_scores = np.zeros(n_reps)

        for u in tqdm(range(0, n_reps)):
            sample_l_n, sample_l_c = lasso_sample(N, C, N1, C1)

            sample, labels = prepare_data(sample_l_n, sample_l_c)

            # prepare the data
            X = sample
            X = X.values
            y = labels
            y = y['Normal'].values


            # save all the predictions
            y_train_all = []
            y_test_all = []
            y_pred_all = []
            y_predictions_all = []
            train_probs_all = []
            train_predictions_all = []


            # cross validation
            skf = RepeatedStratifiedKFold(n_splits=5, n_repeats=10, random_state=50)
            for train_index, test_index in skf.split(X, y):
                X_train, X_test = X[train_index], X[test_index]
                y_train, y_test = y[train_index], y[test_index]

                # scale the data
                scaler = preprocessing.StandardScaler().fit(X_train)
                X_train = scaler.transform(X_train)
                X_test = scaler.transform(X_test)


                # Create the random forest instance
                mod = RandomForestClassifier(n_estimators=20, max_features=0.4, max_depth=10, min_samples_leaf=2)
                # Fit the data
                mod.fit(X_train, y_train)
                # Predict the data
                y_pred = [x[1] for x in mod.predict_proba(X_test)]
                y_predictions = mod.predict(X_test)
                train_predictions = mod.predict(X_train)
                train_probs = [x[1] for x in mod.predict_proba(X_train)]


                # Save all the predictions
                y_pred_all += list(y_pred)
                y_predictions_all += list(y_predictions)
                y_train_all += list(y_train)
                y_test_all += list(y_test)
                train_predictions_all += list(train_predictions)
                train_probs_all += list(train_probs)

            auc_scores[u] = roc_auc_score(y_test_all, y_pred_all)
            accuracy_scores[u] = accuracy_score(y_test_all, y_predictions_all)
            f1_scores[u] = f1_score(y_test_all, y_predictions_all)

        values = pd.DataFrame(columns=['auc',  'accuracy',  'f1'])
        values['auc'] = auc_scores
        values['accuracy'] = accuracy_scores
        values['f1'] = f1_scores
        values.to_csv('./Auc_acc_f1/auc_acc_f1_random_' + str(tumor) + '.csv')

        return values