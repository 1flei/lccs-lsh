#pragma once 

#include "util.h"
#include "pri_queue.h"
#include <string>

//fif::()->unique_ptr<Indexer>
//lff::int top_k->unique_ptr<ListType>
//qf::indexer&->int k->double** query->ListType*->IO
template <class IndexerFactoryFunc, class QueryFunc, class ListFactoryFunc>
void benchmark(
    int qn, 
    const float** query, 
    const Result** results,
    // const std::string& outputFilename,  
    FILE* fp, 
    const IndexerFactoryFunc& fIndexerFactory, 
    const QueryFunc& fQuery, 
    const ListFactoryFunc& fListFactory)
{
    // std::cout << outputFilename << std::endl;
    // FILE* fp = fopen(outputFilename.c_str(), "a+");
    // if (!fp) {
    //     printf("cannot open file %s\n", outputFilename.c_str());
    //     return;
    // }

	timeval start_time, end_time;
    gettimeofday(&start_time, NULL);
    auto lsh = fIndexerFactory();
    gettimeofday(&end_time, NULL);
    float indexing_time = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
    
    printf("Indexing Time: %f Seconds\n\n", indexing_time);
    fprintf(fp, "Indexing Time: %f Seconds\n\n", indexing_time);

    // int kMIPs[] = { 1, 2, 5, 10, 100 };
    int kMIPs[] = { 1, 2, 5, 10};
    int max_round = sizeof(kMIPs) / sizeof(int);
    int top_k = -1;

    float percentage = 0.5;
    // printf("Top-k weighted ANN of AWS_SPHERE_LSH: \n");
    // printf("  Top-k\t\t%.2fRatio\tTime (ms)\tRecall\n", percentage);
    printf("  Top-k\t\tRatio\tTime (ms)\tRecall\n");
    for (int num = 0; num < max_round; num++) {
        gettimeofday(&start_time, NULL);

        top_k = kMIPs[num];
        auto list = fListFactory(top_k);
        std::vector<float> ratios(qn);
        std::vector<float> recalls(qn);

        float overall_ratio = 0.0f;
        float recall = 0.0f;
        for (int i = 0; i < qn; ++i) {
            list->reset();
            fQuery(*lsh, top_k, query[i], list.get());

            recalls[i] = calc_recall(top_k, (const Result*)results[i], list.get());
            recall += calc_recall(top_k, (const Result*)results[i], list.get());

            ratios[i] = calc_ratio(top_k, (const Result*)results[i], list.get());
            overall_ratio += ratios[i];
            // printf("%d, |res|=%d, [%d, %f]; [%d, %f]; %f, %f\n", i, list->size(), list->ith_id(0), list->ith_key(0), results[i][0].id_, results[i][0].key_, recalls[i], ratios[i]);
        }

        gettimeofday(&end_time, NULL);
		float runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec -
			start_time.tv_usec) / 1000000.0f;
        float avg_runtime = (runtime * 1000.0f) / qn;

        sort(ratios.begin(), ratios.end());
        // sort(recalls.begin(), recalls.end());
        float medianRatio = ratios[qn * percentage];
        // float medianRecall = recalls[qn*percentage];

        overall_ratio = overall_ratio / qn;
        recall = recall / qn;

        // printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, medianRatio,
        //     avg_runtime, recall);
        printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio,
            avg_runtime, recall);
        // for (int i = 0; i < qn; i++) {
        //     fprintf(fp, "%d-%f\n", i, ratios[i]);
        // }
        fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, medianRatio, avg_runtime, recall);
    }
    printf("\n");
    fprintf(fp, "\n");
    // fclose(fp);
    return;
}

//fif::()->unique_ptr<Indexer>
//lff::int top_k->unique_ptr<ListType>
//qf::indexer&->int k->double** query->ListType*->IO
template <class IndexerFactoryFunc, class QueryFunc, class ListFactoryFunc>
void benchmark_multiplek(
    int qn, 
    const float** query, 
    const Result** results,
    std::vector<int>& check_ks,
    FILE* fp, 
    const IndexerFactoryFunc& fIndexerFactory, 
    const QueryFunc& fQuery, 
    const ListFactoryFunc& fListFactory)
{
    // std::cout << outputFilename << std::endl;
    // FILE* fp = fopen(outputFilename.c_str(), "a+");
    // if (!fp) {
    //     printf("cannot open file %s\n", outputFilename.c_str());
    //     return;
    // }

	timeval start_time, end_time;
    gettimeofday(&start_time, NULL);
    auto lsh = fIndexerFactory();
    gettimeofday(&end_time, NULL);
    float indexing_time = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec - start_time.tv_usec) / 1000000.0f;
    
    printf("Indexing Time: %f Seconds\n\n", indexing_time);
    fprintf(fp, "Indexing Time: %f Seconds\n\n", indexing_time);

    // printf("memory-usage: %ld (bytes)\n", lsh->get_memory_usage());
    // fprintf(fp, "memory-usage: %ld (bytes)\n", lsh->get_memory_usage());

    // int kMIPs[] = { 1, 2, 5, 10, 100 };
    int kMIPs[] = { 1, 2, 5, 10};
    int max_round = sizeof(kMIPs) / sizeof(int);
    int top_k = -1;


    for(int check_k:check_ks){    
        float percentage = 0.5;

        printf("check_k=%d\n", check_k);
        fprintf(fp, "check_k=%d\n", check_k);

        // printf("Top-k weighted ANN of AWS_SPHERE_LSH: \n");
        // printf("  Top-k\t\t%.2fRatio\tTime (ms)\tRecall\n", percentage);
        printf("  Top-k\t\tRatio\t\tTime (ms)\tRecall\n");
        for (int num = 0; num < max_round; num++) {
            gettimeofday(&start_time, NULL);

            top_k = kMIPs[num];
            auto list = fListFactory(top_k);
            std::vector<float> ratios(qn);
            std::vector<float> recalls(qn);

            float overall_ratio = 0.0f;
            float recall = 0.0f;
            for (int i = 0; i < qn; ++i) {
                list->reset();
                fQuery(*lsh, top_k, check_k, query[i], list.get());

                recalls[i] = calc_recall(top_k, (const Result*)results[i], list.get());
                recall += calc_recall(top_k, (const Result*)results[i], list.get());

                ratios[i] = calc_ratio(top_k, (const Result*)results[i], list.get());
                overall_ratio += ratios[i];
                // printf("%d, |res|=%d, [%d, %f]; [%d, %f]; %f, %f\n", i, list->size(), list->ith_id(0), list->ith_key(0), results[i][0].id_, results[i][0].key_, recalls[i], ratios[i]);
            }

            gettimeofday(&end_time, NULL);
            float runtime = end_time.tv_sec - start_time.tv_sec + (end_time.tv_usec -
                start_time.tv_usec) / 1000000.0f;
            float avg_runtime = (runtime * 1000.0f) / qn;

            sort(ratios.begin(), ratios.end());
            // sort(recalls.begin(), recalls.end());
            float medianRatio = ratios[qn * percentage];
            // float medianRecall = recalls[qn*percentage];

            overall_ratio = overall_ratio / qn;
            recall = recall / qn;

            // printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, medianRatio,
            //     avg_runtime, recall);
            printf("  %3d\t\t%.4f\t\t%.4f\t\t%.2f%%\n", top_k, overall_ratio,
                avg_runtime, recall);
            // for (int i = 0; i < qn; i++) {
            //     fprintf(fp, "%d-%f\n", i, ratios[i]);
            // }
            fprintf(fp, "%d\t%f\t%f\t%f\n", top_k, medianRatio, avg_runtime, recall);
        }
        printf("\n");
        fprintf(fp, "\n");
    }

    // fclose(fp);
    return;
}

template <class IndexerFactoryFunc, class QueryFunc>
void benchmarkMinklist(
    int qn,
    const float** query,
    const Result** results,
    FILE* fp, 
    const IndexerFactoryFunc& fIndexerFactory, 
    const QueryFunc& fQuery)
{
    const auto& fListFactory = [&](int top_k){
        return std::make_unique<MinK_List>(top_k);
    };
    // benchmark(qn, query, results, outputFilename, fIndexerFactory, fQuery, fListFactory);
    benchmark(qn, query, results, fp, fIndexerFactory, fQuery, fListFactory);
}

template <class IndexerFactoryFunc, class QueryFunc>
void benchmarkMaxklist(
    int qn,
    const float** query,
    const Result** results,
    FILE* fp, 
    const IndexerFactoryFunc& fIndexerFactory, 
    const QueryFunc& fQuery)
{
    const auto& fListFactory = [&](int top_k){
        return std::make_unique<MaxK_List>(top_k);
    };
    benchmark(qn, query, results, fp, fIndexerFactory, fQuery, fListFactory);
}

template <class IndexerFactoryFunc, class QueryFunc>
void benchmarkMinklist(
    int qn,
    const float** query,
    const Result** results,
    std::vector<int>& check_ks,
    FILE* fp, 
    const IndexerFactoryFunc& fIndexerFactory, 
    const QueryFunc& fQuery)
{
    const auto& fListFactory = [&](int top_k){
        return std::make_unique<MinK_List>(top_k);
    };
    // benchmark(qn, query, results, outputFilename, fIndexerFactory, fQuery, fListFactory);
    benchmark_multiplek(qn, query, results, check_ks, fp, fIndexerFactory, fQuery, fListFactory);
}

template <class IndexerFactoryFunc, class QueryFunc>
void benchmarkMaxklist(
    int qn,
    const float** query,
    const Result** results,
    std::vector<int>& check_ks,
    FILE* fp, 
    const IndexerFactoryFunc& fIndexerFactory, 
    const QueryFunc& fQuery)
{
    const auto& fListFactory = [&](int top_k){
        return std::make_unique<MaxK_List>(top_k);
    };
    benchmark_multiplek(qn, query, results, check_ks, fp, fIndexerFactory, fQuery, fListFactory);
}