// graph_3color_backtrack.c
// Complete backtracking algorithm for 3-coloring with pre-coloring support
// GUARANTEED to always give the correct answer
// Serves as ground truth for testing other algorithms
// Modified to accept pre-coloring as input and output coloring to a file
//
// Compile: gcc -O3 graph_3color_backtrack.c -o graph_3color_backtrack
// Usage:   ./graph_3color_backtrack input.graph [precoloring.txt]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

typedef struct {
    int N, M;
    int *deg;
    int **nbrs;
} Graph;

// Global statistics
static long long nodes_explored = 0;
static long long backtracks = 0;
static int precolored_count = 0;

// Read graph from file
static Graph *read_graph(const char *path) {
    FILE *f = fopen(path, "r");
    if (!f) {
        perror("fopen");
        return NULL;
    }

    int N, M;
    if (fscanf(f, "%d %d", &N, &M) != 2) {
        fclose(f);
        return NULL;
    }

    Graph *g = calloc(1, sizeof *g);
    g->N = N;
    g->M = M;
    g->deg = calloc(N+1, sizeof(int));

    int *U = malloc(M * sizeof(int));
    int *V = malloc(M * sizeof(int));

    // Read edges and count degrees
    for(int i = 0; i < M; i++) {
        if(fscanf(f, "%d %d", &U[i], &V[i]) != 2) {
            free(U); free(V); free(g->deg); free(g);
            fclose(f);
            return NULL;
        }
        g->deg[U[i]]++;
        g->deg[V[i]]++;
    }

    // Allocate adjacency lists
    g->nbrs = malloc((N+1) * sizeof(int*));
    for(int i = 1; i <= N; i++) {
        g->nbrs[i] = malloc(g->deg[i] * sizeof(int));
        g->deg[i] = 0; // Reset for filling
    }

    // Fill adjacency lists
    for(int i = 0; i < M; i++) {
        int a = U[i], b = V[i];
        g->nbrs[a][g->deg[a]++] = b;
        g->nbrs[b][g->deg[b]++] = a;
    }

    free(U);
    free(V);
    fclose(f);
    return g;
}

// Read pre-coloring from file
static int read_precoloring(const char *path, int *coloring, int N) {
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "Warning: Could not open pre-coloring file %s\n", path);
        return 0;
    }

    char line[256];
    int count = 0;

    while (fgets(line, sizeof(line), f)) {
        // Skip comments and empty lines
        if (line[0] == '#' || line[0] == '\n' || line[0] == '\r') {
            continue;
        }

        int vertex, color;
        if (sscanf(line, "%d %d", &vertex, &color) == 2) {
            if (vertex >= 1 && vertex <= N && color >= 0 && color <= 2) {
                coloring[vertex] = color;
                count++;
            } else {
                fprintf(stderr, "Warning: Invalid pre-coloring entry: vertex %d color %d\n", vertex, color);
            }
        }
    }

    fclose(f);
    printf("Loaded %d pre-colored vertices from %s\n", count, path);
    return count;
}

// Verify pre-coloring is valid
static int verify_precoloring(Graph *g, int *coloring) {
    for(int v = 1; v <= g->N; v++) {
        if(coloring[v] == -1) continue; // Not pre-colored

        if(coloring[v] < 0 || coloring[v] > 2) {
            printf("Invalid pre-color %d for vertex %d\n", coloring[v], v);
            return 0;
        }

        for(int i = 0; i < g->deg[v]; i++) {
            int neighbor = g->nbrs[v][i];
            if(coloring[neighbor] != -1 && coloring[v] == coloring[neighbor]) {
                printf("Pre-coloring conflict: vertices %d and %d both have color %d\n",
                       v, neighbor, coloring[v]);
                return 0;
            }
        }
    }
    return 1;
}

static void free_graph(Graph *g) {
    if(!g) return;
    for(int i = 1; i <= g->N; i++) {
        free(g->nbrs[i]);
    }
    free(g->nbrs);
    free(g->deg);
    free(g);
}

// Check if assigning color c to vertex v is valid
static int is_valid(Graph *g, int *coloring, int v, int c) {
    for(int i = 0; i < g->deg[v]; i++) {
        int neighbor = g->nbrs[v][i];
        if(coloring[neighbor] == c) {
            return 0; // Conflict with neighbor
        }
    }
    return 1;
}

// Recursive backtracking function
static int backtrack(Graph *g, int *coloring, int vertex) {
    nodes_explored++;

    // Base case: all vertices colored
    if(vertex > g->N) {
        return 1; // Success!
    }

    // Skip pre-colored vertices
    if(coloring[vertex] != -1) {
        return backtrack(g, coloring, vertex + 1);
    }

    // Try each color for current vertex
    for(int color = 0; color < 3; color++) {
        if(is_valid(g, coloring, vertex, color)) {
            // Assign color
            coloring[vertex] = color;

            // Recurse to next vertex
            if(backtrack(g, coloring, vertex + 1)) {
                return 1; // Found solution
            }

            // Backtrack
            coloring[vertex] = -1;
            backtracks++;
        }
    }

    return 0; // No valid coloring found
}

// Enhanced backtracking with optimizations
static int enhanced_backtrack(Graph *g, int *coloring, int vertex) {
    nodes_explored++;

    // Base case: all vertices colored
    if(vertex > g->N) {
        return 1;
    }

    // Skip already colored vertices (including pre-colored)
    if(coloring[vertex] != -1) {
        return enhanced_backtrack(g, coloring, vertex + 1);
    }

    // Count available colors for this vertex
    int available[3] = {1, 1, 1}; // Initially all colors available
    int num_available = 3;

    for(int i = 0; i < g->deg[vertex]; i++) {
        int neighbor = g->nbrs[vertex][i];
        int neighbor_color = coloring[neighbor];
        if(neighbor_color != -1) {
            if(available[neighbor_color]) {
                available[neighbor_color] = 0;
                num_available--;
            }
        }
    }

    // Early termination: no colors available
    if(num_available == 0) {
        backtracks++;
        return 0;
    }

    // Try available colors
    for(int color = 0; color < 3; color++) {
        if(available[color]) {
            coloring[vertex] = color;

            if(enhanced_backtrack(g, coloring, vertex + 1)) {
                return 1;
            }

            coloring[vertex] = -1;
            backtracks++;
        }
    }

    return 0;
}

// Choose next uncolored vertex using MRV (Most Constrained Variable) heuristic
static int choose_vertex_mrv(Graph *g, int *coloring) {
    int best_vertex = -1;
    int min_remaining = 4; // More than 3 colors

    for(int v = 1; v <= g->N; v++) {
        if(coloring[v] != -1) continue; // Already colored (including pre-colored)

        // Count remaining colors for vertex v
        int available[3] = {1, 1, 1};
        int remaining = 3;

        for(int i = 0; i < g->deg[v]; i++) {
            int neighbor = g->nbrs[v][i];
            int neighbor_color = coloring[neighbor];
            if(neighbor_color != -1 && available[neighbor_color]) {
                available[neighbor_color] = 0;
                remaining--;
            }
        }

        if(remaining < min_remaining) {
            min_remaining = remaining;
            best_vertex = v;
        }

        // Early termination: found vertex with no valid colors
        if(remaining == 0) {
            return v;
        }
    }

    return best_vertex;
}

// Advanced backtracking with MRV heuristic
static int advanced_backtrack(Graph *g, int *coloring) {
    nodes_explored++;

    // Choose next vertex to color using MRV heuristic
    int vertex = choose_vertex_mrv(g, coloring);

    if(vertex == -1) {
        return 1; // All vertices colored successfully
    }

    // Count available colors
    int available[3] = {1, 1, 1};
    int num_available = 3;

    for(int i = 0; i < g->deg[vertex]; i++) {
        int neighbor = g->nbrs[vertex][i];
        int neighbor_color = coloring[neighbor];
        if(neighbor_color != -1 && available[neighbor_color]) {
            available[neighbor_color] = 0;
            num_available--;
        }
    }

    if(num_available == 0) {
        backtracks++;
        return 0; // No colors available
    }

    // Try each available color
    for(int color = 0; color < 3; color++) {
        if(available[color]) {
            coloring[vertex] = color;

            if(advanced_backtrack(g, coloring)) {
                return 1;
            }

            coloring[vertex] = -1;
            backtracks++;
        }
    }

    return 0;
}

// Verify a complete coloring
static int verify_coloring(Graph *g, int *coloring) {
    for(int v = 1; v <= g->N; v++) {
        if(coloring[v] < 0 || coloring[v] > 2) {
            printf("Invalid color %d for vertex %d\n", coloring[v], v);
            return 0;
        }

        for(int i = 0; i < g->deg[v]; i++) {
            int neighbor = g->nbrs[v][i];
            if(coloring[v] == coloring[neighbor]) {
                printf("Color conflict: vertices %d and %d both have color %d\n",
                       v, neighbor, coloring[v]);
                return 0;
            }
        }
    }
    return 1;
}

// Print coloring (for debugging)
static void print_coloring(Graph *g, int *coloring) {
    printf("Coloring: ");
    for(int v = 1; v <= g->N; v++) {
        printf("%d:%d ", v, coloring[v]);
    }
    printf("\n");
}

// Write coloring to file
static void write_coloring_to_file(const char *graph_filename, Graph *g, int *coloring) {
    // Create output filename by appending .coloring to the input filename
    int len = strlen(graph_filename);
    char *output_filename = malloc(len + 10);
    strcpy(output_filename, graph_filename);
    strcat(output_filename, ".coloring");

    FILE *f = fopen(output_filename, "w");
    if (!f) {
        fprintf(stderr, "Error: Could not create coloring file %s\n", output_filename);
        free(output_filename);
        return;
    }

    // Write header
    fprintf(f, "# 3-Coloring for graph: %s\n", graph_filename);
    fprintf(f, "# Number of vertices: %d\n", g->N);
    if (precolored_count > 0) {
        fprintf(f, "# Pre-colored vertices: %d\n", precolored_count);
    }
    fprintf(f, "# Color mapping: 0=RED, 1=GREEN, 2=BLUE\n");
    fprintf(f, "\n");

    // Write coloring
    for(int v = 1; v <= g->N; v++) {
        fprintf(f, "%d %d\n", v, coloring[v]);
    }

    fclose(f);
    printf("Coloring written to: %s\n", output_filename);

    // Also show color distribution
    int color_count[3] = {0, 0, 0};
    for(int v = 1; v <= g->N; v++) {
        color_count[coloring[v]]++;
    }

    printf("Color distribution:\n");
    const char *color_names[] = {"RED", "GREEN", "BLUE"};
    for(int c = 0; c < 3; c++) {
        double percentage = (color_count[c] * 100.0) / g->N;
        printf("  %s (%d):   %d vertices (%.1f%%)\n",
               color_names[c], c, color_count[c], percentage);
    }

    free(output_filename);
}

int main(int argc, char **argv) {
    if(argc != 2 && argc != 3) {
        fprintf(stderr, "Usage: %s graph.txt [precoloring.txt]\n", argv[0]);
        fprintf(stderr, "  graph.txt      - input graph file\n");
        fprintf(stderr, "  precoloring.txt - optional pre-coloring file\n");
        return 1;
    }

    Graph *g = read_graph(argv[1]);
    if(!g) {
        fprintf(stderr, "Failed to read graph\n");
        return 2;
    }

    printf("Graph: %d vertices, %d edges\n", g->N, g->M);

    // Initialize coloring (-1 = uncolored)
    int *coloring = malloc((g->N + 1) * sizeof(int));
    for(int i = 0; i <= g->N; i++) {
        coloring[i] = -1;
    }

    // Load pre-coloring if provided
    if(argc == 3) {
        precolored_count = read_precoloring(argv[2], coloring, g->N);

        if(precolored_count > 0) {
            // Verify pre-coloring is valid
            if(!verify_precoloring(g, coloring)) {
                fprintf(stderr, "Error: Pre-coloring is invalid (has conflicts)\n");
                free(coloring);
                free_graph(g);
                return 3;
            }
            printf("Pre-coloring verified - no conflicts\n");

            // Show how many vertices remain to be colored
            int remaining = g->N - precolored_count;
            printf("Vertices to color: %d (%.1f%% of total)\n",
                   remaining, (remaining * 100.0) / g->N);
        }
    }

    // Reset statistics
    nodes_explored = 0;
    backtracks = 0;

    clock_t start = clock();

    // Try advanced backtracking first (usually faster)
    int result = advanced_backtrack(g, coloring);

    clock_t end = clock();
    double time_taken = ((double)(end - start)) / CLOCKS_PER_SEC;

    if(result) {
        printf("3-COLORABLE\n");

        // Verify the solution
        if(verify_coloring(g, coloring)) {
            printf("Solution verified!\n");

            // Write coloring to file
            write_coloring_to_file(argv[1], g, coloring);

            // Optionally print the coloring to console (uncomment if needed)
            // print_coloring(g, coloring);
        } else {
            printf("ERROR: Invalid solution found!\n");
        }
    } else {
        printf("UNCOLORED\n");
        printf("The graph cannot be 3-colored");
        if(precolored_count > 0) {
            printf(" with the given pre-coloring");
        }
        printf("\n");
    }

    printf("\nStatistics:\n");
    printf("  Pre-colored vertices: %d\n", precolored_count);
    printf("  Nodes explored: %lld\n", nodes_explored);
    printf("  Backtracks: %lld\n", backtracks);
    printf("  Time taken: %.6f seconds\n", time_taken);

    free(coloring);
    free_graph(g);
    return 0;
}
