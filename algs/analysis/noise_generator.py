from src.models.noise.noise_configuration import NoiseModel


def main():
    #model_cte = NoiseModel("model_constant", "constant", process_noise=1)
    model_rcte = NoiseModel("model_random_constant", "random_constant", process_noise=[1,2,3])
    #model_wn = NoiseModel("model_white_noise", "white_noise", process_noise=1)
    #model_gm = NoiseModel("model_gauss_markov", "gauss_markov", process_noise=1, correlation_time=1)
    #model_rw = NoiseModel("model_random_walk", "random_walk", process_noise=1)


    noise_rcte = model_rcte.gen(100)

    print(noise_rcte)



print("#--------------------------------------------------#")
print("#           Welcome to GNSSNavPy Program           #")
print("#--------------------------------------------------#\n")

if __name__ == "__main__":
    main()
