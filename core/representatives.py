import numpy as np
import numpy.typing as npt
import tqdm
import sys
from sklearn.metrics import silhouette_score
from sklearn.metrics import silhouette_samples
from sklearn.cluster import KMeans

# local
from core.cluster import Cluster
import core.tools as tools


# Implemented by Felipe V. Calderan
def calculate_silhouette(
    data: npt.NDArray[float], labels: npt.NDArray[int], k: int
) -> float:
    """
    Compute Silhouette for each sample, take the average for each cluster and
    return (#clusters that surpassed the general Silhouette) / (#clusters)
    """
    # Calculate silhouette for individual samples
    silh_samples = silhouette_samples(data, labels)
    # Calculate general silhouette
    silh_score = silhouette_score(data, labels)

    # Create cluster-silhouette array and order it by label number
    csilh = np.array(list(zip(labels, silh_samples)))
    csilh = csilh[csilh[:, 0].argsort()]

    # Group silhouette scores belonging in the same cluster together
    groups = np.split(
        csilh[:, 1], np.unique(csilh[:, 0], return_index=True)[1][1:]
    )

    # Return (#clusters that surpassed the general Silhouette) / (#clusters)
    return sum([group.mean() >= silh_score for group in groups]) / k


# Implemented by Felipe V. Calderan
def perform_clustering_n_random(krange: list, 
                                seed_list: int,
                                runs: int, 
                                data: npt.NDArray[float]):
    """
    Perform K-Means clustering n_random times for all possible values of K
    """

    # this stores clustering data of all seeds for all ks
    all_clusterings = []

    # maxk = min(math.sqrt(len(data))+1, krange[1])+1
    maxk = min(len(data)/2, krange[1])+1
    # perform K-Means for all possible values of K
    # for k in range(3, min(math.ceil(len(data)/2), maxsamples)):
    if maxk > len(data):
        maxk = len(data)
    set_k = range(int(krange[0]), int(maxk), int(krange[2]))
    for k in tqdm.tqdm(set_k, desc='Searching K', leave=False, position=1):
        # print("K: ", k, len(seed_list))
        # this stores clustering data for result values of all available seeds
        clustering = []

        # perform K-Means for all available seeds
        for i, seed in enumerate(seed_list):
            kmeans = KMeans(n_clusters=k, n_init=runs, random_state=seed).fit(data)
            score = calculate_silhouette(data, kmeans.labels_, k)
            # print(k, i, score)
            wcss = kmeans.inertia_

            c = Cluster(1)
            c.best_k = k
            c.insert(k, kmeans.labels_, kmeans.cluster_centers_, score, wcss)
            # print("K: ", k, c.labels[0])#kmeans.labels_)
            clustering.append(c)
            # print("Cluster: ",k,i," - ", c, "\n\n")
            

        all_clusterings.append(clustering)
    return all_clusterings, set_k

# Implemented by Felipe V. Calderan
def extract_best_k(krange: list, 
                   clusters_data: list, 
                   set_k: list) -> tuple[int, float, dict]:
    """
    Get the best K available
    """
    k_scores = []

    # compute means and variances
    for i, clus in enumerate(clusters_data):
        scores = []
        for d in clus:
            scores.append(d.scores[0])
        k_scores.append(
            {"k": set_k[i], "mean": np.mean(scores), "var": np.var(scores)}
        )
    # compute the best K
    best_k = [0]
    best_mean = [float("-inf")]

    for s in k_scores:
        # print("KKK: ", s["k"], s["mean"])
        # if there is a better value, keep only it
        if s["mean"] > best_mean[0]:
            best_mean = [s["mean"]]
            best_k = [s["k"]]
        # otherwise, add to the list
        elif s["mean"] == best_mean[0]:
            best_mean.append(s["mean"])
            best_k.append(s["k"])

    # # if there are multiple ks with the same mean, take the largest
    if len(best_k) > 1:
        best_k = int(np.max(best_k))
        best_mean = best_mean[0]
    # othewise, the best k is already the only one available
    else:
        best_k = best_k[0]
        best_mean = best_mean[0]

    # print("Selected: ", best_k)

    return best_k, best_mean, k_scores


# Implemented by Felipe V. Calderan
def pick_best_candidate(krange: list,
                        clusters_data: list[list[Cluster]], 
                        best_k: int) -> Cluster:
    """
    Pick the best candidate from the clusterings that have the best K
    """

    start_k = krange[0]


    Ks = [clus[0].best_k for clus in clusters_data]
    for id, clus in enumerate(clusters_data):
        if clus[0].best_k == best_k:
            idc = id
            best_cluster_data = clus
            break

    best_score = float("-inf")
    best_i = -1

    for i, bcd in enumerate(best_cluster_data):
        if bcd.scores[0] > best_score:
            best_score = bcd.scores[0]
            best_i = i

    return best_cluster_data[best_i]

def get_representatives(params : dict, coulomb, energies, pfiles: str):
    runs = params["RUN_KMEANS"]
    cluster_labels = []
    cluster_centroids = []            

    if len(set(energies)) != 1:
        energies = np.array([n if n is not None else sys.maxsize for n in energies])

    if len(params["KMEANS"])==1:
        k = params["KMEANS"][0]
        print("K   : ", k)
        if k > len(coulomb):
            k = len(coulomb)
        print("KKKK: ", k)
        kmeans = KMeans(n_clusters=k, random_state=0, n_init=runs).fit(coulomb)
        cluster_labels = kmeans.labels_
        cluster_centroids = kmeans.cluster_centers_
    else: 
        # Optimize K with Silhouette
        # print("\t\tSilhouette: k searching")
        pseudo_seeds = [np.random.randint(999999) for i in range(10)]
        evaluation, set_k = perform_clustering_n_random(params["KMEANS"], pseudo_seeds, runs, coulomb)
        if len(evaluation)==0:
            print("\n\n\t\t\tFail - there is not enough samples to perform data clustering")
            try:
                folder_in = params["TMP_FOLDER"]+"/filtered/"
            except:
                folder_in = params["INPUT_FOLDER"]
            print(f"\t\t\tXYZ Files are available at {folder_in}\n\n")
            exit()

        best_k, best_mean, k_scores = extract_best_k(params["KMEANS"], evaluation, set_k)
        # print("Bestk: ", evaluation)
        result = pick_best_candidate(params["KMEANS"], evaluation, best_k)
        # print("\t\tResult from clustering analysis: K = ",best_k)
        cluster_labels = result.labels[0] 
        cluster_centroids = result.centroids[0]

    sel_samples = []
    n_samples = params["MAXSAMPLES"]-1
    file_ids = np.array(range(len(pfiles)))

    for cluster in range(cluster_centroids.shape[0]):
        cluster_points = coulomb[cluster_labels == cluster]
        ids = file_ids[cluster_labels==cluster]
        distances = np.sqrt(np.sum((cluster_points - cluster_centroids[cluster])**2, axis=1))
        min_dist_index = np.argmin(distances)
        closest_id = ids[min_dist_index]
        sel_samples.append(closest_id)
        ids = ids[ids != closest_id]
        if n_samples > 0:
            if n_samples < len(cluster_points)-1:
                samples = ids[np.random.choice(ids.shape[0], 
                                                          n_samples, 
                                                          replace=False)]
            else:
                samples = ids

            sel_samples.extend(samples)

    return sel_samples

