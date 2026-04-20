import { Toaster as ChakraToaster, createToaster, Toast } from "@chakra-ui/react";

export const toaster = createToaster({
  placement: "top",
  duration: 3000,
});

export function Toaster() {
  return (
    <ChakraToaster toaster={toaster}>
      {(toast) => <Toast.Root key={toast.id}><Toast.Title>{toast.title}</Toast.Title><Toast.Description>{toast.description}</Toast.Description></Toast.Root>}
    </ChakraToaster>
  );
}
