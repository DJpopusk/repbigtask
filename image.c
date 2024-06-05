#include "lodepng.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <time.h>

typedef struct Node {
    unsigned char r, g, b, a;
    struct Node *up, *down, *left, *right;
    struct Node *parent;
    int rank;
} Node;

Node* find(Node* x) {
    if (x->parent != x) {
        x->parent = find(x->parent);
    }
    return x->parent;
}
void blur_filter(Node* nodes, int width, int height) {
    int kernel[3][3] = {
        {1, 2, 1},
        {2, 4, 2},
        {1, 2, 1}
    };
    int kernel_sum = 16;

    Node* new_nodes = malloc(width * height * sizeof(Node));
    if (!new_nodes) return;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            int r = 0, g = 0, b = 0;
            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    int ny = y + ky;
                    int nx = x + kx;
                    if (ny >= 0 && ny < height && nx >= 0 && nx < width) {
                        Node* neighbor = &nodes[ny * width + nx];
                        int kvalue = kernel[ky + 1][kx + 1];
                        r += neighbor->r * kvalue;
                        g += neighbor->g * kvalue;
                        b += neighbor->b * kvalue;
                    }
                }
            }
            Node* new_node = &new_nodes[y * width + x];
            new_node->r = r / kernel_sum;
            new_node->g = g / kernel_sum;
            new_node->b = b / kernel_sum;
            new_node->a = 255;

            new_node->up = (y > 0) ? &new_nodes[(y - 1) * width + x] : NULL;
            new_node->down = (y < height - 1) ? &new_nodes[(y + 1) * width + x] : NULL;
            new_node->left = (x > 0) ? &new_nodes[y * width + (x - 1)] : NULL;
            new_node->right = (x < width - 1) ? &new_nodes[y * width + (x + 1)] : NULL;
            new_node->parent = new_node;
            new_node->rank = 0;
        }
    }

    memcpy(nodes, new_nodes, width * height * sizeof(Node));
    free(new_nodes);
}

void sobel_filter(Node* nodes, unsigned w, unsigned h) {
    int x, y;
    int Gx[3][3] = {
        {-1, 0, 1},
        {-2, 0, 2},
        {-1, 0, 1}
    };
    int Gy[3][3] = {
        {1, 2, 1},
        {0, 0, 0},
        {-1, -2, -1}
    };

    Node* tempNodes = (Node*)malloc(w * h * sizeof(Node));

    for (y = 1; y < h - 1; y++) {
        for (x = 1; x < w - 1; x++) {
            float sumX = 0.0;
            float sumY = 0.0;
            Node* current = &nodes[y * w + x];

            for (int i = -1; i <= 1; i++) {
                for (int j = -1; j <= 1; j++) {
                    Node* neighbor = &nodes[(y + i) * w + (x + j)];
                    unsigned char r = neighbor->r;
                    sumX += Gx[i + 1][j + 1] * r;
                    sumY += Gy[i + 1][j + 1] * r;
                }
            }

            int sum = (int)sqrt(sumX * sumX + sumY * sumY);
            if (sum > 255) sum = 255;
            if (sum < 0) sum = 0;

            tempNodes[y * w + x].r = sum;
            tempNodes[y * w + x].g = sum;
            tempNodes[y * w + x].b = sum;
            tempNodes[y * w + x].a = 255;
        }
    }

    for (y = 0; y < h; y++) {
        for (x = 0; x < w; x++) {
            nodes[y * w + x].r = tempNodes[y * w + x].r;
            nodes[y * w + x].g = tempNodes[y * w + x].g;
            nodes[y * w + x].b = tempNodes[y * w + x].b;
            nodes[y * w + x].a = tempNodes[y * w + x].a;
        }
    }

    free(tempNodes);
}


int compare(const void *a, const void *b) {
    return (*(unsigned char*)a - *(unsigned char*)b);
}

void median_filter(Node* nodes, int width, int height) {
    Node* new_nodes = malloc(width * height * sizeof(Node));
    if (!new_nodes) return;

    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            unsigned char r_values[9], g_values[9], b_values[9];
            int k = 0;
            for (int ky = -1; ky <= 1; ky++) {
                for (int kx = -1; kx <= 1; kx++) {
                    int ny = y + ky;
                    int nx = x + kx;
                    if (ny >= 0 && ny < height && nx >= 0 && nx < width) {
                        Node* neighbor = &nodes[ny * width + nx];
                        r_values[k] = neighbor->r;
                        g_values[k] = neighbor->g;
                        b_values[k] = neighbor->b;
                        k++;
                    }
                }
            }
            qsort(r_values, k, sizeof(unsigned char), compare);
            qsort(g_values, k, sizeof(unsigned char), compare);
            qsort(b_values, k, sizeof(unsigned char), compare);

            Node* new_node = &new_nodes[y * width + x];
            new_node->r = r_values[k / 2];
            new_node->g = g_values[k / 2];
            new_node->b = b_values[k / 2];
            new_node->a = 255;

            new_node->up = (y > 0) ? &new_nodes[(y - 1) * width + x] : NULL;
            new_node->down = (y < height - 1) ? &new_nodes[(y + 1) * width + x] : NULL;
            new_node->left = (x > 0) ? &new_nodes[y * width + (x - 1)] : NULL;
            new_node->right = (x < width - 1) ? &new_nodes[y * width + (x + 1)] : NULL;
            new_node->parent = new_node;
            new_node->rank = 0;
        }
    }

    memcpy(nodes, new_nodes, width * height * sizeof(Node));
    free(new_nodes);
}


void union_set(Node* x, Node* y, double epsilon) {
    if (x->r < 60) {
        return;
    }
    Node* px = find(x);
    Node* py = find(y);

    double color_difference = sqrt(pow(x->r - y->r, 2) + pow(x->g - y->g, 2) + pow(x->b - y->b, 2));
    if (px != py && color_difference < epsilon) {
        if (px->rank > py->rank) {
            py->parent = px;
        } else {
            px->parent = py;
            if (px->rank == py->rank) {
                py->rank++;
            }
        }
    }
}

Node* create_graph(const char *filename, int *width, int *height) {
    unsigned char *image = NULL;
    int error = lodepng_decode32_file(&image, width, height, filename);
    if (error) {
        printf("error %u: %s\n", error, lodepng_error_text(error));
        return NULL;
    }
  for (int i = 0; i < *height * *width * 4; i += 4) {
    char gray = (image[i + 0] +image[i + 1] + image[i + 2]) / 3;
    image[i + 0] = gray;
    image[i + 1] = gray;
    image[i + 2] = gray;
    image[i + 3] = image[i + 3];
  }

    Node* nodes = malloc(*width * *height * sizeof(Node));
    if (!nodes) {
        free(image);
        return NULL;
    }
    
    for (unsigned y = 0; y < *height; ++y) {
        for (unsigned x = 0; x < *width; ++x) {
            Node* node = &nodes[y * *width + x];
            unsigned char* pixel = &image[(y * *width + x) * 4];
            node->r = pixel[0];
            node->g = pixel[1];
            node->b = pixel[2];
            node->a = pixel[3];
            node->up = y > 0 ? &nodes[(y - 1) * *width + x] : NULL;
            node->down = y < *height - 1 ? &nodes[(y + 1) * *width + x] : NULL;
            node->left = x > 0 ? &nodes[y * *width + (x - 1)] : NULL;
            node->right = x < *width - 1 ? &nodes[y * *width + (x + 1)] : NULL;
            node->parent = node;
            node->rank = 0;
        }
    }


    free(image);
    return nodes;
}

void find_components(Node* nodes, int width, int height, double epsilon) {
    for (int y = 0; y < height; y++) {
        for (int x = 0; x < width; x++) {
            Node* node = &nodes[y * width + x];
            if (node->up) {
                union_set(node, node->up, epsilon);
            }
            if (node->down) {
                union_set(node, node->down, epsilon);
            }
            if (node->left) {
                union_set(node, node->left, epsilon);
            }
            if (node->right) {
                union_set(node, node->right, epsilon);
            }
        }
    }
}

void free_graph(Node* nodes) {
    free(nodes);
}

void color_components_and_count(Node* nodes, int width, int height) {
    unsigned char* output_image = malloc(width * height * 4 * sizeof(unsigned char));
    int* component_sizes = calloc(width * height, sizeof(int));
    int total_components = 0;

    srand(time(NULL));
    for (int i = 0; i < width * height; i++) {
        Node* p = find(&nodes[i]);
        if (p == &nodes[i]) {
            if (component_sizes[i] < 1) {
                p->r = 0;
                p->g = 0;
                p->b = 0;
            } else {
                p->r = rand() % 256;
                p->g = rand() % 256;
                p->b = rand() % 256;
            }
            total_components++;
        }
        output_image[4 * i + 0] = p->r;
        output_image[4 * i + 1] = p->g;
        output_image[4 * i + 2] = p->b;
        output_image[4 * i + 3] = 255;
        component_sizes[p - nodes]++;
    }

    char *output_filename = "output1.png";
    lodepng_encode32_file(output_filename, output_image, width, height);


    free(output_image);
    free(component_sizes);
}
int main(int argc, char* argv[]) {
    if (argc != 2) {
        fprintf(stderr, "Usage: %s <filename>\n", argv[0]);
        return 1;
    }

    int width, height;
    char *filename = argv[1];
    Node* nodes = create_graph(filename, &width, &height);
    if (!nodes) {
        fprintf(stderr, "Error creating graph from file: %s\n", filename);
        return 1;
    }
    double epsilon = 15.0;
    find_components(nodes, width, height, epsilon);
    color_components_and_count(nodes, width, height);
    blur_filter(nodes, width, height);
    median_filter(nodes, width, height);
    sobel_filter(nodes, width, height);
    free_graph(nodes);
    return 0;
}
