//very important line
localStorage.removeItem("animationPlayed"); // comment this line if you dont want to loop
//VVIP line

if (!localStorage.getItem("animationPlayed")) {

  const animation = lottie.loadAnimation({
    container: document.getElementById("animation"),
    renderer: "svg",
    loop: false,//play only once for false, multiple for true
    autoplay: true,
    path: "https://dulran.com/Logo-5-remixver2.json"// path to the json file for the animation
  });

  animation.addEventListener("complete", () => {
    localStorage.setItem("animationPlayed", "true");// this should be true if I want the loading screen to show up once for the same computer but
    //but I set it equal to flase for now so i can see the logo showing up
    
    const animDiv = document.getElementById("animation");
    animDiv.style.transition = "opacity 0.8s ease";
    animDiv.style.opacity = "0";
  
    setTimeout(() => {
      animDiv.remove();
    }, 800);
 });//hiding the animation and removing the layer after one loop

} else {
  document.getElementById("animation").style.display = "none";
}
