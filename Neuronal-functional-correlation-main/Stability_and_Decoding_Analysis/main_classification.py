from scipy.io import savemat
import h5py
import sys
sys.path.insert(0, '//166.111.72.183/slfm/WLF/cebra_test_try_GPU//third_party')
def noliner_by_trails_RAW_with_label_for_every_trials_RUSH(file_path,sti_kind,param, trace_kind,classify_kind,corr_win):
    import random
    import numpy as np
    import torch
    from sklearn.neighbors import KNeighborsClassifier
    import matplotlib
    matplotlib.use('TkAgg')

    if torch.cuda.is_available():
        DEVICE = "cuda"
    else:
        DEVICE = "cpu"

    # ################################################################################################################## get matrix fun
    def get_multi_neurons_trials_matrix_by_sti(trace, stimuli, param):
        before = int(param['show_sti_before'] * param['fs'])
        during = int(param['sti_during_time'] * param['fs'])
        after = int(param['show_sti_after'] * param['fs'])
        start_edge = stimuli['start_edge']
        end_edge = stimuli['end_edge']
        sti_labels = np.unique(stimuli['stimuili_label_ind'], return_index=False)
        sti_num = len(sti_labels)
        # sti_num = len(sti_labels) - 1
        trials_cell = {}
        # for j in range(sti_num):
        # for j in range(1,sti_num):
        for j in range(1,sti_num):
            sti_i = int(sti_labels[j])
            trials_matrix_for_each_sti = {}
            trials_num = 0
            for i in range(len(start_edge)):
                try:
                    if stimuli['stimuili_label_ind'][i] == sti_i:
                        trial = trace[int(start_edge[i] - before): int(start_edge[i] + during + after), :]
                        trials_matrix_for_each_sti[trials_num] = np.array(trial)
                        trials_num = trials_num + 1
                except:
                    if stimuli['stimuili_label_ind'][:, i] == sti_i:
                        trial = trace[int(start_edge[i] - before): int(start_edge[i] + during + after), :]
                        trials_matrix_for_each_sti[trials_num] = np.array(trial)
                        trials_num = trials_num + 1
            trials_cell[sti_i] = trials_matrix_for_each_sti
        return trials_cell, before, during, after, sti_num, trials_num

    def nolinear_for_trials(trace,stimuli_mat,sti_kind,test_index, classify_kind):
        neuron_num = trace.shape[1]
        # ################################################################################################################## get matrix
        [trials_cell, before, during, after, sti_num, trials_num] = get_multi_neurons_trials_matrix_by_sti(trace[:,:], stimuli, param)
        # trials_num = trials_num - 1
        # ################################################################################################################## set train & test
        train_trials = {}
        train_trials_label = []
        test_trials = {}
        test_trials_label = []

        train_trials_index = [i for i in range(trials_num) if i != test_index]
        test_trials_index = test_index
        train_trials_index = train_trials_index[7:] # pre 8 trials used for select neurons

        train_trials_i = 0
        # for i in range(sti_num):
        for i in range(1,sti_num):
            for j in train_trials_index: # change here
                train_trials[train_trials_i] = trials_cell[i+1][j]
                train_trials_label.append(i+1)
                train_trials_i = train_trials_i+1
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> get trial matrix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ')

        test_trials_i = 0
        # for i in range(sti_num):
        for i in range(1,sti_num):
            test_trials[test_trials_i] = trials_cell[i+1][test_trials_index]
            test_trials_label.append(i+1)
            test_trials_i = test_trials_i+1
        print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> get trial matrix >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> ')

        # ################################################################################################################## shuffle train & test
        randnum = random.randint(0,100)
        random.seed(randnum)
        random.shuffle(train_trials_label)
        random.seed(randnum)
        random.shuffle(train_trials)
        print('shuffled train data')

        randnum = random.randint(0,100)
        random.seed(randnum)
        random.shuffle(test_trials_label)
        random.seed(randnum)
        random.shuffle(test_trials)
        print('shuffled test data')

        # ################################################################################################################## classify method
        def knn_decode(train_trials, train_trials_label, test_trials, test_trials_label):
            train_data = []
            for key, value in train_trials.items():
                train_data.append(value)
            train_data = np.array(train_data)
            train_data = train_data.reshape(train_data.shape[0], -1)
            train_labels = np.array(train_trials_label)

            test_data = []
            for key, value in test_trials.items():
                test_data.append(value)
            test_data = np.array(test_data)
            test_data = test_data.reshape(test_data.shape[0], -1)
            test_labels = np.array(test_trials_label)

            params = np.power(np.linspace(1, 5, 5, dtype=int), 2)
            # params = np.power(np.linspace(1, 3, 3, dtype=int), 2)
            errs = []
            for n in params:
                train_decoder = KNeighborsClassifier(n_neighbors=n, metric='cosine')
                train_valid_idx = int(len(train_labels) / 9 * 8)
                train_decoder.fit(train_data[:train_valid_idx], train_labels[:train_valid_idx])
                pred = train_decoder.predict(train_data[train_valid_idx:])
                err = train_labels[train_valid_idx:] - pred
                errs.append(abs(err).sum())

            train_decoder = KNeighborsClassifier(n_neighbors=params[np.argmin(errs)], metric='cosine')
            train_decoder.fit(train_data, train_labels)

            accuracy, test_label, pred = train_decoder.score(test_data, test_labels)
            # print(test_label.T)
            # print(pred)
            # print("Accuracy:", accuracy)
            # print(f'RAW my data all neurons: {accuracy*100:.2f}%')
            return accuracy, test_label, pred

        from sklearn.svm import SVC
        def svm_decode(train_trials, train_trials_label, test_trials, test_trials_label):
            train_data = []
            for _, value in train_trials.items():
                train_data.append(value)
            train_data = np.array(train_data)
            train_data = train_data.reshape(train_data.shape[0], -1)
            train_labels = np.array(train_trials_label)

            test_data = []
            for _, value in test_trials.items():
                test_data.append(value)
            test_data = np.array(test_data)
            test_data = test_data.reshape(test_data.shape[0], -1)
            test_labels = np.array(test_trials_label)

            params = np.logspace(-3, 3, 7)
            errs = []
            best_score = 0
            best_c = 1
            for c in params:
                svc_classifier = SVC(C=c, kernel='rbf')
                train_valid_idx = int(len(train_labels) / 9 * 8)
                svc_classifier.fit(train_data[:train_valid_idx], train_labels[:train_valid_idx])
                score = svc_classifier.score(train_data[train_valid_idx:], train_labels[train_valid_idx:])
                score = score[0]

                errs.append(1 - score)
                if score > best_score:
                    best_score = score
                    best_c = c

            svc_classifier = SVC(C=best_c, kernel='rbf')
            svc_classifier.fit(train_data, train_labels)

            accuracy = svc_classifier.score(test_data, test_labels)
            pred = svc_classifier.predict(test_data)

            return accuracy, test_labels, pred

        from sklearn.ensemble import BaggingClassifier
        def ensemble_svm_decode(train_trials, train_trials_label, test_trials, test_trials_label,neuron_num):
            train_data = np.array([value for _, value in train_trials.items()])
            train_data = train_data.reshape(train_data.shape[0], -1)
            train_labels = np.array(train_trials_label)

            test_data = np.array([value for _, value in test_trials.items()])
            test_data = test_data.reshape(test_data.shape[0], -1)
            test_labels = np.array(test_trials_label)

            # base_clf = SVC(C=1)
            # # base_clf = SVC(C=c, kernel='rbf')

            params = np.logspace(-3, 3, 7)
            errs = []
            best_score = 0
            best_c = 1

            for c in params:
                svc_classifier = SVC(C=c, kernel='rbf')
                train_valid_idx = int(len(train_labels) / 9 * 8)
                svc_classifier.fit(train_data[:train_valid_idx], train_labels[:train_valid_idx])
                score = svc_classifier.score(train_data[train_valid_idx:], train_labels[train_valid_idx:])
                score = score[0]

                errs.append(1 - score)
                if score > best_score:
                    best_score = score
                    best_c = c

            base_clf = SVC(C=best_c, kernel='rbf')
            n_estimators = 300
            # max_features = 'auto'

            bagging_clf = BaggingClassifier(base_clf, n_estimators, max_samples=0.8, max_features=0.8)
            bagging_clf.fit(train_data, train_labels)

            accuracy = bagging_clf.score(test_data, test_labels)
            preds = bagging_clf.predict(test_data)

            return accuracy, test_labels, preds

        from sklearn.ensemble import RandomForestClassifier
        def rf_decode(train_trials, train_trials_label, test_trials, test_trials_label,neuron_num):
            # 将训练数据集整理为合适的形式
            train_data = []
            for _, value in train_trials.items():
                train_data.append(value)
            train_data = np.array(train_data)
            train_data = train_data.reshape(train_data.shape[0], -1)
            train_labels = np.array(train_trials_label)

            # 将测试数据集整理为合适的形式
            test_data = []
            for _, value in test_trials.items():
                test_data.append(value)
            test_data = np.array(test_data)
            test_data = test_data.reshape(test_data.shape[0], -1)
            test_labels = np.array(test_trials_label)

            # 设置随机森林的参数，可以根据需要调整这些参数
            n_estimators = neuron_num*1 # 树的数量
            max_features = 'auto'  # 寻找最佳分割时要考虑的特征数量（'auto' 代表所有特征）

            train_decoder = RandomForestClassifier(n_estimators=n_estimators, max_features=max_features, random_state=0)
            train_decoder.fit(train_data, train_labels)

            # 评估模型性能
            accuracy = train_decoder.score(test_data, test_labels)
            pred = train_decoder.predict(test_data)

            return accuracy, test_labels, pred

        if classify_kind == 'rf':
            acc, test_label, pred = rf_decode(train_trials, train_trials_label, test_trials, test_trials_label,neuron_num)
            acc = acc[0]

        if classify_kind == 'knn':
            acc, test_label, pred = knn_decode(train_trials, train_trials_label, test_trials, test_trials_label)

        if classify_kind == 'svm':
            acc, test_label, pred = svm_decode(train_trials, train_trials_label, test_trials, test_trials_label)
            acc = acc[0]

        if classify_kind == 'ensemble_svm':
            acc, test_label, pred = ensemble_svm_decode(train_trials, train_trials_label, test_trials, test_trials_label,neuron_num)
            acc = acc[0]


        print(f'RAW my data all neurons: {acc*100:.2f}%')
        print(test_label)
        print(pred)
        print('contrast')
        return acc, test_label, pred

    # ################################################################################################################## load data
    if corr_win == '-2':
        data_mat = h5py.File(file_path + '/top_250_edges_trial_variance.mat', 'r')
    if corr_win == '0':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel.mat', 'r')
    if corr_win == '2':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win2.mat', 'r')
    if corr_win == '4':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win4.mat', 'r')
    if corr_win == '8':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win8.mat', 'r')
    if corr_win == '12':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win12.mat', 'r')
    if corr_win == '16':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win16.mat', 'r')
    if corr_win == '20':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win20.mat', 'r')
    if corr_win == '24':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win24.mat', 'r')
    if corr_win == '28':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win28.mat', 'r')
    if corr_win == '32':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win32.mat', 'r')
    if corr_win == '36':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win36.mat', 'r')
    if corr_win == '40':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win40.mat', 'r')
    if corr_win == '44':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win44.mat', 'r')
    if corr_win == '48':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win48.mat', 'r')
    if corr_win == '52':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win52.mat', 'r')
    if corr_win == '56':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win56.mat', 'r')
    if corr_win == '60':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win60.mat', 'r')
    if corr_win == '64':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win64.mat', 'r')
    if corr_win == '68':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win68.mat', 'r')
    if corr_win == '72':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win72.mat', 'r')
    if corr_win == '76':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win76.mat', 'r')
    if corr_win == '80':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win80.mat', 'r')
    if corr_win == '96':
        data_mat = h5py.File(file_path + '/neuron_region_lr_with_select_neurons_withnozerolabel_win96.mat', 'r')

    # if corr_win == 'top1%':
    #     data_mat = h5py.File(file_path + '/top_0.01_edges_trial_variance.mat', 'r')
    # if corr_win == 'top5%':
    #     data_mat = h5py.File(file_path + '/top_0.05_edges_trial_variance.mat', 'r')
    # if corr_win == 'top10%':
    #     data_mat = h5py.File(file_path + '/top_0.10_edges_trial_variance.mat', 'r')
    # if corr_win == 'top20%':
    #     data_mat = h5py.File(file_path + '/top_0.20_edges_trial_variance.mat', 'r')
    # if corr_win == 'top30%':
    #     data_mat = h5py.File(file_path + '/top_0.30_edges_trial_variance.mat', 'r')
    # if corr_win == 'top40%':
    #     data_mat = h5py.File(file_path + '/top_0.40_edges_trial_variance.mat', 'r')
    # if corr_win == 'top50%':
    #     data_mat = h5py.File(file_path + '/top_0.50_edges_trial_variance.mat', 'r')
    # if corr_win == 'top60%':
    #     data_mat = h5py.File(file_path + '/top_0.60_edges_trial_variance.mat', 'r')
    # if corr_win == 'top70%':
    #     data_mat = h5py.File(file_path + '/top_0.70_edges_trial_variance.mat', 'r')
    # if corr_win == 'top80%':
    #     data_mat = h5py.File(file_path + '/top_0.80_edges_trial_variance.mat', 'r')
    # if corr_win == 'top90%':
    #     data_mat = h5py.File(file_path + '/top_0.90_edges_trial_variance.mat', 'r')
    # if corr_win == 'top100%':
    #     data_mat = h5py.File(file_path + '/top_1.00_edges_trial_variance.mat', 'r')

    if corr_win == 'top0.1%':
        data_mat = h5py.File(file_path + '/top_0.001_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.2%':
        data_mat = h5py.File(file_path + '/top_0.002_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.3%':
        data_mat = h5py.File(file_path + '/top_0.003_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.4%':
        data_mat = h5py.File(file_path + '/top_0.004_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.5%':
        data_mat = h5py.File(file_path + '/top_0.005_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.6%':
        data_mat = h5py.File(file_path + '/top_0.006_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.7%':
        data_mat = h5py.File(file_path + '/top_0.007_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.8%':
        data_mat = h5py.File(file_path + '/top_0.008_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.9%':
        data_mat = h5py.File(file_path + '/top_0.009_edges_trial_variance.mat', 'r')
    if corr_win == 'top1.5%':
        data_mat = h5py.File(file_path + '/top_0.015_edges_trial_variance.mat', 'r')
    if corr_win == 'top1%':
        data_mat = h5py.File(file_path + '/top_0.010_edges_trial_variance.mat', 'r')
    if corr_win == 'top2%':
        data_mat = h5py.File(file_path + '/top_0.020_edges_trial_variance.mat', 'r')
    if corr_win == 'top3%':
        data_mat = h5py.File(file_path + '/top_0.030_edges_trial_variance.mat', 'r')
    if corr_win == 'top4%':
        data_mat = h5py.File(file_path + '/top_0.040_edges_trial_variance.mat', 'r')
    if corr_win == 'top5%':
        data_mat = h5py.File(file_path + '/top_0.050_edges_trial_variance.mat', 'r')
    if corr_win == 'top6%':
        data_mat = h5py.File(file_path + '/top_0.060_edges_trial_variance.mat', 'r')
    if corr_win == 'top7%':
        data_mat = h5py.File(file_path + '/top_0.070_edges_trial_variance.mat', 'r')
    if corr_win == 'top8%':
        data_mat = h5py.File(file_path + '/top_0.080_edges_trial_variance.mat', 'r')
    if corr_win == 'top9%':
        data_mat = h5py.File(file_path + '/top_0.090_edges_trial_variance.mat', 'r')
    if corr_win == 'top10%':
        data_mat = h5py.File(file_path + '/top_0.100_edges_trial_variance.mat', 'r')
    if corr_win == 'top20%':
        data_mat = h5py.File(file_path + '/top_0.200_edges_trial_variance.mat', 'r')
    if corr_win == 'top30%':
        data_mat = h5py.File(file_path + '/top_0.300_edges_trial_variance.mat', 'r')
    if corr_win == 'top40%':
        data_mat = h5py.File(file_path + '/top_0.400_edges_trial_variance.mat', 'r')
    if corr_win == 'top50%':
        data_mat = h5py.File(file_path + '/top_0.500_edges_trial_variance.mat', 'r')
    if corr_win == 'top60%':
        data_mat = h5py.File(file_path + '/top_0.600_edges_trial_variance.mat', 'r')
    if corr_win == 'top70%':
        data_mat = h5py.File(file_path + '/top_0.700_edges_trial_variance.mat', 'r')
    if corr_win == 'top80%':
        data_mat = h5py.File(file_path + '/top_0.800_edges_trial_variance.mat', 'r')
    if corr_win == 'top90%':
        data_mat = h5py.File(file_path + '/top_0.900_edges_trial_variance.mat', 'r')
    if corr_win == 'top100%':
        data_mat = h5py.File(file_path + '/top_1.000_edges_trial_variance.mat', 'r')

    if corr_win == 'top0.010%':
        data_mat = h5py.File(file_path + '/top_0.0001_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.020%':
        data_mat = h5py.File(file_path + '/top_0.0002_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.040%':
        data_mat = h5py.File(file_path + '/top_0.0004_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.080%':
        data_mat = h5py.File(file_path + '/top_0.0008_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.160%':
        data_mat = h5py.File(file_path + '/top_0.0016_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.320%':
        data_mat = h5py.File(file_path + '/top_0.0032_edges_trial_variance.mat', 'r')
    if corr_win == 'top0.640%':
        data_mat = h5py.File(file_path + '/top_0.0064_edges_trial_variance.mat', 'r')
    if corr_win == 'top1.280%':
        data_mat = h5py.File(file_path + '/top_0.0128_edges_trial_variance.mat', 'r')
    if corr_win == 'top2.560%':
        data_mat = h5py.File(file_path + '/top_0.0256_edges_trial_variance.mat', 'r')
    if corr_win == 'top5.120%':
        data_mat = h5py.File(file_path + '/top_0.0512_edges_trial_variance.mat', 'r')
    if corr_win == 'top10.240%':
        data_mat = h5py.File(file_path + '/top_0.1024_edges_trial_variance.mat', 'r')
    if corr_win == 'top20.480%':
        data_mat = h5py.File(file_path + '/top_0.2048_edges_trial_variance.mat', 'r')
    if corr_win == 'top40.960%':
        data_mat = h5py.File(file_path + '/top_0.4096_edges_trial_variance.mat', 'r')
    if corr_win == 'top81.920%':
        data_mat = h5py.File(file_path + '/top_0.8192_edges_trial_variance.mat', 'r')
    if corr_win == 'top100.000%':
        data_mat = h5py.File(file_path + '/top_1.0000_edges_trial_variance.mat', 'r')

    if corr_win == 'topnum10':
        data_mat = h5py.File(file_path + '/top_10_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum20':
        data_mat = h5py.File(file_path + '/top_20_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum30':
        data_mat = h5py.File(file_path + '/top_30_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum50':
        data_mat = h5py.File(file_path + '/top_50_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum100':
        data_mat = h5py.File(file_path + '/top_100_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum150':
        data_mat = h5py.File(file_path + '/top_150_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum200':
        data_mat = h5py.File(file_path + '/top_200_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum250':
        data_mat = h5py.File(file_path + '/top_250_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum300':
        data_mat = h5py.File(file_path + '/top_300_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum400':
        data_mat = h5py.File(file_path + '/top_400_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum500':
        data_mat = h5py.File(file_path + '/top_500_edges_trial_variance.mat', 'r')
    if corr_win == 'topnum600':
        data_mat = h5py.File(file_path + '/top_600_edges_trial_variance.mat', 'r')

    if corr_win == 'topTV_edge_group1':
        data_mat = h5py.File(file_path + '/edge_stage1.mat', 'r')
    if corr_win == 'topTV_edge_group2':
        data_mat = h5py.File(file_path + '/edge_stage2.mat', 'r')
    if corr_win == 'topTV_edge_group3':
        data_mat = h5py.File(file_path + '/edge_stage3.mat', 'r')
    if corr_win == 'topTV_edge_group4':
        data_mat = h5py.File(file_path + '/edge_stage4.mat', 'r')
    if corr_win == 'topTV_edge_group5':
        data_mat = h5py.File(file_path + '/edge_stage5.mat', 'r')
    if corr_win == 'topTV_edge_group6':
        data_mat = h5py.File(file_path + '/edge_stage6.mat', 'r')
    if corr_win == 'topTV_edge_group7':
        data_mat = h5py.File(file_path + '/edge_stage7.mat', 'r')
    if corr_win == 'topTV_edge_group8':
        data_mat = h5py.File(file_path + '/edge_stage8.mat', 'r')
    if corr_win == 'topTV_edge_group9':
        data_mat = h5py.File(file_path + '/edge_stage9.mat', 'r')
    if corr_win == 'topTV_edge_group10':
        data_mat = h5py.File(file_path + '/edge_stage10.mat', 'r')
    if corr_win == 'shuffle_TV_edge_group':
        data_mat = h5py.File(file_path + '/edge_shuffle.mat', 'r')
    if corr_win == 'topTV_neuron_group1':
        data_mat = h5py.File(file_path + '/all_neuron_stage1.mat', 'r')
    if corr_win == 'topTV_neuron_group2':
        data_mat = h5py.File(file_path + '/all_neuron_stage2.mat', 'r')
    if corr_win == 'topTV_neuron_group3':
        data_mat = h5py.File(file_path + '/all_neuron_stage3.mat', 'r')
    if corr_win == 'topTV_neuron_group4':
        data_mat = h5py.File(file_path + '/all_neuron_stage4.mat', 'r')
    if corr_win == 'topTV_neuron_group5':
        data_mat = h5py.File(file_path + '/all_neuron_stage5.mat', 'r')
    if corr_win == 'topTV_neuron_group6':
        data_mat = h5py.File(file_path + '/all_neuron_stage6.mat', 'r')
    if corr_win == 'topTV_neuron_group7':
        data_mat = h5py.File(file_path + '/all_neuron_stage7.mat', 'r')
    if corr_win == 'topTV_neuron_group8':
        data_mat = h5py.File(file_path + '/all_neuron_stage8.mat', 'r')
    if corr_win == 'topTV_neuron_group9':
        data_mat = h5py.File(file_path + '/all_neuron_stage9.mat', 'r')
    if corr_win == 'topTV_neuron_group10':
        data_mat = h5py.File(file_path + '/all_neuron_stage10.mat', 'r')
    if corr_win == 'shuffle_TV_neuron_group':
        data_mat = h5py.File(file_path + '/all_neuron_shuffle.mat', 'r')

    if corr_win == 'threshold_neuron_group1':
        data_mat = h5py.File(file_path + '/all_neuron_threshold1.mat', 'r')
    if corr_win == 'threshold_neuron_group2':
        data_mat = h5py.File(file_path + '/all_neuron_threshold2.mat', 'r')
    if corr_win == 'threshold_neuron_group3':
        data_mat = h5py.File(file_path + '/all_neuron_threshold3.mat', 'r')
    if corr_win == 'threshold_neuron_group4':
        data_mat = h5py.File(file_path + '/all_neuron_threshold4.mat', 'r')
    if corr_win == 'threshold_neuron_group5':
        data_mat = h5py.File(file_path + '/all_neuron_threshold5.mat', 'r')
    if corr_win == 'threshold_neuron_group6':
        data_mat = h5py.File(file_path + '/all_neuron_threshold6.mat', 'r')
    if corr_win == 'threshold_neuron_group7':
        data_mat = h5py.File(file_path + '/all_neuron_threshold7.mat', 'r')
    if corr_win == 'threshold_neuron_group8':
        data_mat = h5py.File(file_path + '/all_neuron_threshold8.mat', 'r')

    if corr_win == 'threshold_edge_group1':
        data_mat = h5py.File(file_path + '/all_edge_threshold1.mat', 'r')
    if corr_win == 'threshold_edge_group2':
        data_mat = h5py.File(file_path + '/all_edge_threshold2.mat', 'r')
    if corr_win == 'threshold_edge_group3':
        data_mat = h5py.File(file_path + '/all_edge_threshold3.mat', 'r')
    if corr_win == 'threshold_edge_group4':
        data_mat = h5py.File(file_path + '/all_edge_threshold4.mat', 'r')
    if corr_win == 'threshold_edge_group5':
        data_mat = h5py.File(file_path + '/all_edge_threshold5.mat', 'r')
    if corr_win == 'threshold_edge_group6':
        data_mat = h5py.File(file_path + '/all_edge_threshold6.mat', 'r')
    if corr_win == 'threshold_edge_group7':
        data_mat = h5py.File(file_path + '/all_edge_threshold7.mat', 'r')
    if corr_win == 'threshold_edge_group8':
        data_mat = h5py.File(file_path + '/all_edge_threshold8.mat', 'r')

    if corr_win == 'higher_threshold_edge_group1':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold1.mat', 'r')
    if corr_win == 'higher_threshold_edge_group2':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold2.mat', 'r')
    if corr_win == 'higher_threshold_edge_group3':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold3.mat', 'r')
    if corr_win == 'higher_threshold_edge_group4':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold4.mat', 'r')
    if corr_win == 'higher_threshold_edge_group5':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold5.mat', 'r')
    if corr_win == 'higher_threshold_edge_group6':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold6.mat', 'r')
    if corr_win == 'higher_threshold_edge_group7':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold7.mat', 'r')
    if corr_win == 'higher_threshold_edge_group8':
        data_mat = h5py.File(file_path + '/all_edge_act_deact_higher_threshold8.mat', 'r')

    data = {key: value for key, value in data_mat.items() if not key.startswith('__')}
    print('loaded cal data')

    # ################################################################################################################## load stimuli
    stimuli_mat = h5py.File(file_path+'/visual_stimuli_with_label.mat', 'r')
    stimuli = {key: value for key, value in stimuli_mat.items() if not key.startswith('__')}
    print('load stimuli data')

    trace = data[trace_kind]
    neuron_num = trace.shape[1]

    all_acc = []
    all_label = []
    all_pred = []
    for test_index in range(8, 20):
    # for test_index in range(7, 19):
    # for test_index in range(5, 9):
        RAW_acc, RAW_label, RAW_pred = nolinear_for_trials(trace, stimuli_mat, sti_kind, test_index, classify_kind)
        all_acc.append(RAW_acc)
        all_label.append(RAW_label)
        all_pred.append(RAW_pred)
    mean_acc = np.mean(all_acc)
    print(f'>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> mean_acc: {mean_acc * 100:.2f}%')

    return mean_acc,all_acc,all_label,all_pred,neuron_num

if __name__ == '__main__':
    file_path = '//166.111.72.183/slfm/WLF/20250101/1_20250101_M742_AI162_natureimages_IPD3sISI10s_5image20trial_20X_C2/result_for_thre8' #ok

    # ##########################################################  set param ############################################
    print('>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>> loaded parm >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>')
    readme = 'AI:all neuron intensity; SIa:selected neuron intensity active;  ' \
             'SCa:selected neuron correlation active;  ' \
             'SCKa:key selected neuron correlation active; ' \
             'SIad:selected neuron intensity active & deactive;  ' \
             'SCad:selected neuron correlation active & deactive;  ' \
             'SCKad:key selected neuron correlation active & deactive'
    sti_kind = 'natureimage'

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  choose method
    param = {'fs': 4,
             'show_sti_before': 3,  # change here according to different sti
             'sti_during_time': 3,  # change here according to different sti
             'show_sti_after': 3,  # change here according to different sti
             'delay_after': 0}
    param_for_corr = {'fs': 4,
             'show_sti_before': 3,  # change here according to different sti
             'sti_during_time': 3,  # change here according to different sti
             'show_sti_after': 0,  # change here according to different sti
             'delay_after': 0}

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  choose method
    # classify_kind = 'knn'
    classify_kind = 'svm'
    # classify_kind = 'rf'

    # >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>  base classification example
    AI_mean_acc, AI_all_acc, AI_all_label, AI_all_pred, AI_neuron_num = \
        noliner_by_trails_RAW_with_label_for_every_trials_RUSH(file_path, sti_kind, param, 'valid_C', classify_kind,'0')
    SIa_mean_acc, SIa_all_acc, SIa_all_label, SIa_all_pred, SIa_neuron_num = \
        noliner_by_trails_RAW_with_label_for_every_trials_RUSH(file_path, sti_kind, param, 'select_neurons_C', classify_kind,'0')
    SCa_mean_acc, SCa_all_acc, SCa_all_label, SCa_all_pred, SCa_neuron_num = \
        noliner_by_trails_RAW_with_label_for_every_trials_RUSH(file_path, sti_kind, param_for_corr, 'select_neurons_C_sort_by_region_corr_matrix', classify_kind,'12')

    savemat(file_path + '/' + sti_kind + classify_kind + '_classify_result.mat',
            {'AI_mean_acc': AI_mean_acc, 'AI_all_acc': AI_all_acc, 'AI_all_label': AI_all_label,'AI_all_pred': AI_all_pred, 'AI_neuron_num': AI_neuron_num,
             'SIa_mean_acc': SIa_mean_acc, 'SIa_all_acc': SIa_all_acc, 'SIa_all_label': SIa_all_label,'SIa_all_pred': SIa_all_pred, 'SIa_neuron_num': SIa_neuron_num,
             'SCa_mean_acc': SCa_mean_acc, 'SCa_all_acc': SCa_all_acc, 'SCa_all_label': SCa_all_label,'SCa_all_pred': SCa_all_pred, 'SCa_neuron_num': SCa_neuron_num, 'readme': readme})

