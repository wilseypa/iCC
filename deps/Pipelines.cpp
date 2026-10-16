#include "Pipelines.hpp"

#include <iostream>
#include <utility>

#include "DimensionFrame.hpp"
#include "SimplexUtility.hpp"

namespace
{
void printSimplexLabels(
    const BinomialTable& binomial_table,
    const char* const prefix,
    const SimplexBindex simplex_bindex,
    const std::size_t total_label_count,
    const std::size_t simplex_dimension)
{
    std::cout << "    " << prefix << " labels: ";
    if (simplex_bindex < 0)
    {
        std::cout << "(none)\n";
        return;
    }

    const auto labels = SimplexUtility::getSimplexVertices(
        binomial_table,
        simplex_bindex,
        total_label_count,
        simplex_dimension);

    if (labels.empty())
    {
        std::cout << "(empty)";
    }
    else
    {
        for (std::size_t i = 0; i < labels.size(); ++i)
        {
            if (i > 0)
                std::cout << ' ';
            std::cout << labels[i];
        }
    }
    std::cout << '\n';
}

void printBirthDeath(const PersistentPairInfo& pair)
{
    std::cout << pair.facet_weight << ", " << pair.cofacet_weight;
}

void printPersistencePair(
    const PersistentPairInfo& pair,
    const BinomialTable& binomial_table,
    const std::size_t dimension,
    const std::size_t total_label_count)
{
    std::cout << "  (";
    printBirthDeath(pair);
    std::cout << ")\n";

    printSimplexLabels(
        binomial_table,
        "birth facet",
        pair.facet_bindex,
        total_label_count,
        dimension - 1);
    printSimplexLabels(
        binomial_table,
        "death cofacet",
        pair.cofacet_bindex,
        total_label_count,
        dimension);
}

void printTopSimplex(
    const char* const name,
    const SimplexBindex bindex,
    const std::size_t dimension,
    const PipelineRuntime& runtime,
    const WindowState& window,
    const SimplexEnumerator& simplex_enumerator)
{
    std::cout << "    " << name << " simplex:";
    if (bindex < 0)
    {
        std::cout << " (none)\n";
        return;
    }

    const auto labels = SimplexUtility::getSimplexVertices(
        runtime.binomialTable(), bindex, window.totalLabelCount(), dimension);
    const auto realization = simplex_enumerator.getGeometricSimplexRealization(
        labels, window.pseudoVertices(), window.bounds().eps_hi);

    std::cout << "\n      labels:";
    for (const auto label : labels)
        std::cout << ' ' << label;
    std::cout << "\n      geometric realization (original vertices):";
    for (const auto vertex : realization.vertices)
        std::cout << ' ' << vertex;
    std::cout << '\n';
}

void printPvInformation(
    const WindowState& window,
    const std::size_t new_pv_count)
{
    std::cout << "  new PVs: " << new_pv_count << '\n'
              << "  total PVs: " << window.pseudoVertices().size() << '\n';
    if (window.pseudoVertices().empty())
        std::cout << "  (empty)\n";

    for (std::size_t pv_index = 0;
         pv_index < window.pseudoVertices().size();
         ++pv_index)
    {
        const auto& pv = window.pseudoVertices()[pv_index];
        std::cout << "  PV label " << window.originalVertexCount() + pv_index
                  << ": size=" << pv.representatives.size()
                  << " diameter=" << pv.diameter
                  << "\n    original vertices:";
        for (const auto vertex : pv.representatives)
            std::cout << ' ' << vertex;
        std::cout << '\n';
    }
    std::cout << std::endl;
}
}

PersistentHomologyPipeline::PersistentHomologyPipeline(
    DistanceMatrix distance_matrix,
    PipelineConfig config)
    : runtime_(
          std::move(distance_matrix),
          std::move(config),
          PipelineMode::RegVRPH)
{
}

void PersistentHomologyPipeline::run()
{
    WindowState window(runtime_.distanceMatrix().vertexCount());

    const WindowBounds bounds{
        .eps_lo = 0.0,
        .eps_hi = runtime_.config().finalEpsilon(),
        .is_final = true};
    auto edges = window.prepareWindow(runtime_, bounds);

    std::cout << "total point number: "
              << runtime_.distanceMatrix().vertexCount() << std::endl;

    DimensionFrame frame(runtime_, window, std::move(edges));
    for (;;)
    {
        frame.matchPersistence(false);

        if (frame.dimension() >= 2)
        {
            std::cout << "Processing dim " << frame.dimension()
                      << ", cofacet number: " << frame.cofacetCount()
                      << ", facet number: " << frame.facetCount()
                      << '\n'
                      << "dimensional persistent pairs:" << std::endl;

            if (frame.persistentPairs().empty())
            {
                std::cout << "  (empty)" << std::endl;
            }
            else
            {
                for (const auto& pair : frame.persistentPairs())
                {
                    if (runtime_.config().verbose)
                    {
                        printPersistencePair(
                            pair,
                            runtime_.binomialTable(),
                            frame.dimension(),
                            window.totalLabelCount());
                    }
                    else
                    {
                        std::cout << "  (";
                        printBirthDeath(pair);
                        std::cout << ")\n";
                    }
                }
            }
        }

        if (!frame.advance())
            break;
    }
}

PwphPipeline::PwphPipeline(
    DistanceMatrix distance_matrix,
    PipelineConfig config)
    : runtime_(
          std::move(distance_matrix),
          std::move(config),
          PipelineMode::PwPH),
      window_state_(runtime_.distanceMatrix().vertexCount())
{
}

void PwphPipeline::run()
{
    const auto& config = runtime_.config();
    for (std::size_t window_index = 0;
         window_index < config.eps_breaks.size();
         ++window_index)
    {
        const WindowBounds bounds{
            .eps_lo = window_index == 0
                ? 0.0
                : config.eps_breaks[window_index - 1],
            .eps_hi = config.eps_breaks[window_index],
            .is_final = window_index + 1 == config.eps_breaks.size()};

        if (!config.verbose)
        {
            std::cout << "window lower bound, window upper bound\n";
            std::cout << bounds.eps_lo << ", " << bounds.eps_hi << '\n';
            std::cout << "dimension, birth weight, death weight\n";
        }
        else
        {
            std::cout << "epsilon window: " << bounds.eps_lo << ' '
                      << bounds.eps_hi << '\n';
        }

        auto edges = window_state_.prepareWindow(runtime_, bounds);
        const SimplexEnumerator simplex_enumerator(
            runtime_.distanceMatrix(), runtime_.binomialTable());

        DependencySupportBatch support_batch;
        {
            DimensionFrame frame(runtime_, window_state_, std::move(edges));
            for (;;)
            {
                if (config.verbose && frame.dimension() >= 2)
                {
                    std::cout << "interface dimension: " << frame.dimension()
                              << "\n  cofacet list size: " << frame.cofacetCount()
                              << "\n  facet list size: " << frame.facetCount()
                              << "\n  active facet list size: " << frame.activeFacetCount()
                              << std::endl;
                }
                else if (!config.verbose && frame.dimension() >= 2 &&
                         frame.isTopDimension())
                {
                    std::cout << "top interface dimension: " << frame.dimension()
                              << "\n  cofacet list size: " << frame.cofacetCount()
                              << "\n  facet list size: " << frame.facetCount()
                              << std::endl;
                }

                frame.matchPersistence(!bounds.is_final);

                if (frame.dimension() >= 2)
                {
                    if (config.verbose)
                    {
                        std::cout << "PH dimension " << frame.dimension() - 1
                                  << ":" << std::endl;
                    }

                    bool printed_any = false;
                    for (const auto& pair : frame.persistentPairs())
                    {
                        if (!(pair.cofacet_weight >= bounds.eps_lo ||
                              pair.cofacet_weight < 0.0))
                        {
                            continue;
                        }

                        if (config.verbose)
                        {
                            std::cout << "  (";
                            printBirthDeath(pair);
                            std::cout << ")\n";
                            if (frame.isTopDimension())
                            {
                                printTopSimplex(
                                    "birth", pair.facet_bindex,
                                    frame.dimension() - 1, runtime_, window_state_,
                                    simplex_enumerator);
                                printTopSimplex(
                                    "death", pair.cofacet_bindex,
                                    frame.dimension(), runtime_, window_state_,
                                    simplex_enumerator);
                            }
                        }
                        else
                        {
                            std::cout << frame.dimension() << ", ";
                            printBirthDeath(pair);
                            std::cout << '\n';
                        }

                        printed_any = true;
                    }

                    if (config.verbose && !printed_any)
                    {
                        std::cout
                            << "  (no new interval or surviving interval from previous eps range)"
                            << std::endl;
                    }
                }

                if (frame.isTopDimension() && !bounds.is_final)
                    support_batch = frame.takeDependencySupportBatch();

                if (!frame.advance())
                    break;
            }
        }

        // DimensionFrame was the final consumer of this window's geometry.
        // Support postprocessing needs persistent labels/PVs and eps_hi only.
        window_state_.invalidateCurrentWindow();

        if (!config.verbose)
            std::cout << '\n';

        if (bounds.is_final)
            break;

        if (config.verbose)
            std::cout << "PV information:\n";

        const std::size_t previous_pv_count =
            window_state_.pseudoVertices().size();
        support_processor_.processPwph(
            runtime_, window_state_, bounds.eps_hi, std::move(support_batch));
        runtime_.ensureBinomialCapacity(window_state_.totalLabelCount());
        const std::size_t new_pv_count =
            window_state_.pseudoVertices().size() - previous_pv_count;

        if (config.verbose)
            printPvInformation(window_state_, new_pv_count);
    }
}

ApparentPairPipeline::ApparentPairPipeline(
    DistanceMatrix distance_matrix,
    PipelineConfig config)
    : runtime_(
          std::move(distance_matrix),
          std::move(config),
          PipelineMode::RegVRPH)
{
}

void ApparentPairPipeline::run()
{
    WindowState window(runtime_.distanceMatrix().vertexCount());

    const WindowBounds bounds{
        .eps_lo = 0.0,
        .eps_hi = runtime_.config().finalEpsilon(),
        .is_final = true};
    auto edges = window.prepareWindow(runtime_, bounds);

    std::cout << "dimension, unmatched weight\n";

    DimensionFrame frame(runtime_, window, std::move(edges));
    for (;;)
    {
        frame.matchApparentPairsOnly();

        for (const double weight : frame.unmatchedFacetWeights())
            std::cout << frame.facetDimension() << ", " << weight << '\n';

        if (frame.isTopDimension() && frame.dimension() >= 2)
        {
            for (const double weight : frame.unmatchedTopCofacetWeights())
                std::cout << frame.dimension() << ", " << weight << '\n';
        }

        if (!frame.advance())
            break;
    }
}
